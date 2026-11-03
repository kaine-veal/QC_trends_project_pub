"""
plot_variant_trends.py

This script loads archived monthly variant QC CSV files and creates visualisations.
This follows the same pattern as plot_trends.py for interop metrics.

Usage workflow:
1. Archive monthly data using archive_variant_qc_month.py
2. Use this script to plot trends from archived data

Key features:
1. Lists available archived months
2. Loads data from one or more monthly archives
3. Creates boxplots with outlier handling
4. Adds visual separators between months when plotting multiple months
5. Validates metrics and provides helpful error messages
"""

import argparse
import re
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

# Logger setup
from utils.logger import create_logger
from utils.config import OUTPUT_ROOT
logger = create_logger('variant_QC.log')

# Archive location for monthly CSVs
ARCHIVE_BASE = OUTPUT_ROOT / "qc_archive_variantqc"

# Output location for plots
PLOT_OUTPUT_DIR = OUTPUT_ROOT / "month_2_output"

# Metric-specific visualisation limits and outlier thresholds
METRIC_LIMITS = {
    "HetHom_ratio": {
        "y_max": 3.5,
        "outlier_threshold": 3.5,
        "description": "Het/Hom ratio (typical range: 1.5-2.5)"
    },
    "InsDel_ratio": {
        "y_max": 3.0,  # True Ins/Del ratio can range from 0 to 3+ (0.2-2.0 typical)
        "outlier_threshold_low": 0.1,  # Flag if < 0.1 (10x more deletions than insertions)
        "outlier_threshold": 3.0,  # Flag if > 3.0 (3x more insertions than deletions)
        "description": "Ins/Del ratio: insertions/deletions (typical range: 0.5-2.0)"
    },
    "TsTv_ratio": {
        "y_max": 5.0,
        "outlier_threshold_low": 1.5,
        "outlier_threshold": 5.0,
        "description": "Ts/Tv ratio (flag values <1.5 or >5.0)"
    },
    "TotalVariants": {
        "y_max": 500,             
        "outlier_threshold_low": 200,
        "outlier_threshold": 500,
        "description": "Total variants (flag <200 or >500)"
    }
}


# ============================================================================
# Argument Parsing
# ============================================================================
def parse_args():
    """
    Parse command line arguments for the script.

    Returns:
        argparse.Namespace: Parsed arguments containing:
            - months: Comma-separated list of archived months to load (YYYY-MM format)
            - metric: Metric to plot
            - out: Output file path for the plot
            - list_months: Flag to list available archived months
    """
    ap = argparse.ArgumentParser(
        description="Plot Panel variant QC metrics from archived data",
    )
    
    # Main options
    ap.add_argument("--months", help="Comma-separated list of archived months in YYYY-MM format (e.g., 2025-01,2025-02)")
    ap.add_argument("--metric", help="Metric to plot (e.g., HetHom_ratio, TsTv_ratio)")
    ap.add_argument("--out", help="Output plot filename (e.g., hethom_trend.png)")
    ap.add_argument("--list-months", action="store_true", help="List all available archived months")
    
    args = ap.parse_args()
    
    # Validate that either --list-months or required arguments are provided
    if not args.list_months and not args.months:
        ap.error("Either --months or --list-months must be provided")
    
    # Validate that --metric and --out are provided when not listing
    if not args.list_months and (not args.metric or not args.out):
        ap.error("--metric and --out are required (unless using --list-months)")
    
    return args


# ============================================================================
# Archive Management Functions
# ============================================================================
def list_available_months():
    """
    List all available archived months in the variant QC archive.
    Displays months organized by year with status indicators.
    """
    if not ARCHIVE_BASE.exists():
        print(f"No archive found at: {ARCHIVE_BASE}")
        print("Run archive_variant_qc_month.py to create monthly archives first.")
        return
    
    # Collect all month directories (must match YYYY-MM format)
    months = sorted([d.name for d in ARCHIVE_BASE.iterdir() 
                     if d.is_dir() and re.match(r'^\d{4}-\d{2}$', d.name)])
    
    if not months:
        print(f"No archived months found in: {ARCHIVE_BASE}")
        return
    
    print(f"\nAvailable archived months (Panel variant QC):")
    print("=" * 60)
    
    # Group by year for organized display
    by_year = {}
    for month in months:
        year = month.split('-')[0]
        if year not in by_year:
            by_year[year] = []
        
        # Check if CSV exists
        csv_path = ARCHIVE_BASE / month / f"variant_qc_panel_{month}.csv"
        status = "✓" if csv_path.exists() else "✗ (missing CSV)"
        by_year[year].append((month, status))
    
    # Display by year
    for year in sorted(by_year.keys()):
        print(f"\n{year}:")
        for month, status in by_year[year]:
            print(f"  {status} {month}")
    
    print("\n" + "=" * 60)
    print(f"Total: {len(months)} archived months")
    print()


def load_archived_month(month: str) -> pd.DataFrame:
    """
    Load variant QC data from an archived month.
    
    Parameters:
        month (str): Month in YYYY-MM format (e.g., "2025-01")
    
    Returns:
        pd.DataFrame: Loaded metrics data with an additional 'archive_month' column
    
    Raises:
        FileNotFoundError: If the archived CSV doesn't exist
    """
    month_dir = ARCHIVE_BASE / month
    csv_path = month_dir / f"variant_qc_panel_{month}.csv"
    
    if not csv_path.exists():
        logger.error(f"Archived CSV not found: {csv_path}")
        raise FileNotFoundError(f"Archived CSV not found: {csv_path}")
    
    logger.info(f"Loading archived data from: {csv_path}")
    df = pd.read_csv(csv_path)
    
    # Add archive_month column to track which month this data came from
    df['archive_month'] = month
    
    return df


def combine_monthly_archives(months: list) -> pd.DataFrame:
    """
    Load and combine data from multiple archived months.
    
    Parameters:
        months (list): List of months in YYYY-MM format
    
    Returns:
        pd.DataFrame: Combined metrics data from all months
    
    Raises:
        SystemExit: If no data can be loaded from any of the specified months
    """
    all_dfs = []
    
    for month in months:
        try:
            df = load_archived_month(month)
            all_dfs.append(df)
            logger.info(f"Loaded {len(df)} rows from {month}")
        except Exception as e:
            logger.error(f"Failed to load month {month}: {e}")
            print(f"Warning: Could not load month {month}. Skipping...")
            continue
    
    if not all_dfs:
        logger.error("No data found for any of the specified months")
        raise SystemExit("Error: No data found for any of the specified months")
    
    # Combine all dataframes
    combined_df = pd.concat(all_dfs, ignore_index=True)
    logger.info(f"Combined {len(months)} months into {len(combined_df)} total rows")
    
    return combined_df


# ============================================================================
# Plotting Function
# ============================================================================
def plot_variant_metric(df, metric, months, output_path):
    """
    Create boxplot visualization for a variant QC metric with outlier handling.
    
    Parameters:
        df (pd.DataFrame): DataFrame with metrics and BTCH columns
        metric (str): Metric to plot
        months (list): List of months being plotted
        output_path (Path): Path to save the plot
    """
    logger.info(f"Creating plot for metric: {metric}")
    
    # Clean data - convert to numeric and drop NaNs
    df[metric] = pd.to_numeric(df[metric], errors="coerce")
    initial_count = len(df)
    df = df.dropna(subset=[metric])
    dropped_count = initial_count - len(df)
    
    if dropped_count > 0:
        logger.warning(f"Dropped {dropped_count} NaN values for {metric}")
    
    if df.empty:
        logger.error(f"No valid data to plot for {metric}")
        raise SystemExit(f"Error: No valid data for metric {metric}")
    
    # Check for metric-specific limits
    has_limits = metric in METRIC_LIMITS
    
    # Initialize variables
    y_max = None
    y_min = 0
    high = None
    low = None
    num_outliers = 0
    df_normal = df
    df_outliers = pd.DataFrame()
    
    if has_limits:
        limits = METRIC_LIMITS[metric]
        y_max = limits.get("y_max")
        high = limits.get("outlier_threshold")
        low = limits.get("outlier_threshold_low")
        
        # Flag outliers
        df["is_outlier"] = False
        if high is not None:
            df.loc[df[metric] > high, "is_outlier"] = True
        if low is not None:
            df.loc[df[metric] < low, "is_outlier"] = True
        
        num_outliers = df["is_outlier"].sum()
        df_normal = df[~df["is_outlier"]]
        df_outliers = df[df["is_outlier"]]
        
        # Set y_min based on low threshold if it exists
        if low is not None:
            y_min = low
    
    # Set up Seaborn style
    sns.set(style="whitegrid", context="notebook")
    plt.figure(figsize=(14, 6))
    
    # Create boxplot (using normal data only if limits exist)
    ax = sns.boxplot(
        data=df_normal,
        x="BTCH",
        y=metric,
        color="lightgray",
        fliersize=0  # Hide outlier dots (we add our own)
    )
    
    # Add vertical separators between runs
    for i in range(len(df["BTCH"].unique())):
        plt.axvline(i - 0.5, color="gray", linestyle="--", linewidth=0.8)
    
    # Overlay dots for normal data
    sns.stripplot(
        data=df_normal,
        x="BTCH",
        y=metric,
        color="black",
        size=3,
        jitter=True,
        alpha=0.7
    )
    
    # Apply y-axis limits if specified
    if has_limits and y_max is not None:
        ax.set_ylim(y_min, y_max)
    
    # Build title
    title = f"{metric} per run (Panel)"
    if len(months) > 1:
        if len(months) <= 3:
            month_str = ", ".join(months)
            title += f" ({month_str})"
        else:
            title += f" ({len(months)} months)"
    
    ax.set_title(title, fontsize=14)
    ax.set_ylabel(metric, fontsize=11)
    ax.set_xlabel("Run (BTCH ID)", fontsize=11)
    plt.xticks(rotation=45, ha="right", fontsize=9)
    
    plt.tight_layout()
    
    # Add outlier annotation if present
    if has_limits and num_outliers > 0:
        outlier_min = df_outliers[metric].min()
        outlier_max = df_outliers[metric].max()
        
        # Build threshold text
        threshold_parts = []
        if low is not None:
            threshold_parts.append(f"<{low}")
        if high is not None:
            threshold_parts.append(f">{high}")
        
        threshold_text = " or ".join(threshold_parts) if threshold_parts else "outside range"
        
        # Add footer
        outlier_info = (
            f"{num_outliers} outlier(s) ({threshold_text}) "
            f"(range: {outlier_min:.2f}-{outlier_max:.2f}) excluded from scale"
        )
        plt.figtext(0.98, 0.02, outlier_info, 
                    ha="right", va="bottom", fontsize=8, style="italic", color="gray")
    
    # Save plot
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Plot successfully saved: {output_path}")
    
    # Print summary
    print(f"\nMetric: {metric}")
    print(f"Total samples: {initial_count}")
    print(f"Valid samples (after NaN removal): {len(df)}")
    print(f"Normal samples: {len(df_normal)}")
    
    if has_limits and num_outliers > 0:
        threshold_desc = []
        if low is not None:
            threshold_desc.append(f"<{low}")
        if high is not None:
            threshold_desc.append(f">{high}")
        threshold_str = " or ".join(threshold_desc)
        
        print(f"Outliers ({threshold_str}): {num_outliers}")
        print(f"Outlier range: {df_outliers[metric].min():.2f} - {df_outliers[metric].max():.2f}")


# ============================================================================
# Main Function
# ============================================================================
def main():
    """
    Main orchestration function that:
    1. Parses command line arguments
    2. Handles --list-months flag
    3. Loads data from archived months
    4. Creates and saves the plot with selected metric
    """
    logger.info("=" * 60)
    logger.info("Starting plot_variant_trends.py")
    logger.info("=" * 60)

    args = parse_args()
    
    # Handle --list-months flag (early exit)
    if args.list_months:
        list_available_months()
        return

    # Log inputs
    logger.info(f"Metric requested: {args.metric}")
    logger.info(f"Output file: {args.out}")
    logger.info("=" * 60)

    # Parse and validate month list
    months = [m.strip() for m in args.months.split(",")]
    
    # Validate month format
    for month in months:
        try:
            datetime.strptime(month, "%Y-%m")
        except ValueError:
            logger.error(f"Invalid month format: {month}")
            raise SystemExit(f"Error: Month '{month}' is not in YYYY-MM format")
    
    # Load data from archived months
    logger.info(f"Loading data from archived months: {', '.join(months)}")
    df = combine_monthly_archives(months)
    
    # Validate BTCH column exists
    if "BTCH" not in df.columns:
        logger.error("CSV missing required 'BTCH' column")
        raise ValueError("CSV must contain a 'BTCH' column for grouping runs.")
    
    # Validate that requested metric exists
    if args.metric not in df.columns:
        available = ', '.join(df.columns)
        logger.error(f"Metric '{args.metric}' not found in data")
        raise ValueError(f"Metric '{args.metric}' not found. Available columns: {available}")

    # Set up output path
    output_path = PLOT_OUTPUT_DIR / args.out
    
    # Create plot
    plot_variant_metric(df, args.metric, months, output_path)

    logger.info("=" * 60)
    logger.info("plot_variant_trends.py completed successfully")
    logger.info("=" * 60)


# ============================================================================
# Entry Point
# ============================================================================
if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.critical(f"Script failed with exception: {e}", exc_info=True)
        raise
