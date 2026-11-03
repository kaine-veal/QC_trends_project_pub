"""
plot_trends.py

This script loads archived monthly CSV files and creates visualisations of Illumina 
sequencing metrics over time. This is the second step in the QC workflow - after 
archiving monthly data with archive_month.py, use this script to plot trends.

Key features:
1. Lists available archived months
2. Loads data from one or more monthly archives
3. Filters by run type (RACP, Genotyping, Panel)
4. Creates scatter plots with threshold lines where appropriate
5. Adds visual separators between months when plotting multiple months
6. Validates metric/runtype compatibility to prevent confusing plots

The script enforces important restrictions:
- Q30_% and PhiX_Error_Rate_% can only be plotted for a SINGLE run type
  (because each run type has different quality thresholds)
- Total_Runtime_hours and Total_Yield_Gb can be plotted for any combination of run types
"""

import re
import argparse
from pathlib import Path
import pandas as pd
import sys
import matplotlib.pyplot as plt
from datetime import datetime


# Logger setup
from utils.logger import create_logger
from utils.config import OUTPUT_ROOT
logger = create_logger('run_QC.log')

# Archive location for monthly CSVs
ARCHIVE_BASE = OUTPUT_ROOT / "qc_archive_runqc"

# Metric aliases for user convenience (users can type "q30" instead of "Q30_%")
METRIC_ALIASES = {
    "q30": "Q30_%",
    "error_rate": "PhiX_Error_Rate_%",
    "phix_error": "PhiX_Error_Rate_%",
    "runtime": "Total_Runtime_hours",
    "runtime_hours": "Total_Runtime_hours",
    "yield": "Total_Yield_Gb"
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
            - metrics: List of metrics to plot (Q30_%, PhiX_Error_Rate_%, etc.)
            - runtype: Optional list of run types to filter (RACP, Genotyping, Panel)
            - out: Output file path for the plot
            - list_months: Flag to list available archived months
    """
    ap = argparse.ArgumentParser(
        description="Plot Illumina sequencing run metrics from archived data",
    )
    # Main options
    ap.add_argument("--months", help="Comma-separated list of archived months in YYYY-MM format (e.g., 2025-01,2025-02)")
    ap.add_argument("--metrics", nargs="+", help="Metrics to plot (e.g., Q30_%% PhiX_Error_Rate_%%)")
    ap.add_argument("--runtype", nargs="*", choices=["RACP", "Genotyping", "Panel"],
                    help="Filter by run type(s)")
    ap.add_argument("--out", help="Output plot file name (will be saved to plot_output/)")
    ap.add_argument("--list-months", action="store_true", help="List all available archived months")
    
    args = ap.parse_args()
    
    # Validate that either --list-months or required arguments are provided
    if not args.list_months and not args.months:
        ap.error("Either --months or --list-months must be provided")
    
    # Validate that --metrics and --out are provided when not listing
    if not args.list_months and (not args.metrics or not args.out):
        ap.error("--metrics and --out are required (unless using --list-months)")
    
    return args


# ============================================================================
# Archive Management Functions
# ============================================================================
def list_available_months():
    """
    List all available archived months in the interop archive.
    Displays months organized by year with status indicators.
    
    For each month, shows:
    - ✓ if CSV exists
    - ✗ if month directory exists but CSV is missing
    """
    if not ARCHIVE_BASE.exists():
        print(f"No archive found at: {ARCHIVE_BASE}")
        print("Run archive_month.py to create monthly archives first.")
        return
    
    # Collect all month directories (must match YYYY-MM format)
    months = sorted([d.name for d in ARCHIVE_BASE.iterdir() 
                     if d.is_dir() and re.match(r'^\d{4}-\d{2}$', d.name)])
    
    if not months:
        print(f"No archived months found in: {ARCHIVE_BASE}")
        return
    
    print(f"\nAvailable archived months:")
    print("=" * 60)
    
    # Group by year for organized display
    by_year = {}
    for month in months:
        year = month.split('-')[0]
        if year not in by_year:
            by_year[year] = []
        
        # Check if CSV exists
        csv_path = ARCHIVE_BASE / month / f"interop_metrics_{month}.csv"
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
    Load metrics data from an archived month.
    
    Parameters:
        month (str): Month in YYYY-MM format (e.g., "2025-01")
    
    Returns:
        pd.DataFrame: Loaded metrics data with an additional 'archive_month' column
                      to track which month this data came from
    
    Raises:
        FileNotFoundError: If the archived CSV doesn't exist
    """
    month_dir = ARCHIVE_BASE / month
    csv_path = month_dir / f"interop_metrics_{month}.csv"
    
    if not csv_path.exists():
        logger.error(f"Archived CSV not found: {csv_path}")
        raise FileNotFoundError(f"Archived CSV not found: {csv_path}")
    
    logger.info(f"Loading archived data from: {csv_path}")
    df = pd.read_csv(csv_path)
    
    # Add archive_month column to track which month this data came from
    # This is used later for adding visual month separators in plots
    df['archive_month'] = month
    
    return df


def combine_monthly_archives(months: list) -> pd.DataFrame:
    """
    Load and combine data from multiple archived months.
    
    This allows plotting trends across multiple months, with visual separators
    between each month on the plot.
    
    Parameters:
        months (list): List of months in YYYY-MM format (e.g., ['2025-01', '2025-02'])
    
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
    
    # Combine all dataframes vertically (one after another)
    combined_df = pd.concat(all_dfs, ignore_index=True)
    logger.info(f"Combined {len(months)} months into {len(combined_df)} total rows")
    
    return combined_df


# ============================================================================
# Metric Validation
# ============================================================================
def validate_metric_compatibility(metrics, runtypes):
    """
    Validate that selected metrics are compatible with the number of runtype selections.
    
    CRITICAL RULE: Some metrics have different thresholds for different run types,
    so they can only be plotted when a single run type is selected.
    
    Runtype-specific metrics (require exactly ONE runtype):
    - Q30_% (different thresholds: Panel=80%, RACP=75%, Genotyping=75%)
    - PhiX_Error_Rate_% (quality standards differ by runtype)
    
    Universal metrics (can be plotted for any combination of run types):
    - Total_Runtime_hours (runtime is comparable across all run types)
    - Total_Yield_Gb (yield is comparable across all run types)
    
    Parameters:
        metrics (list): Selected metrics to plot.
        runtypes (list): Selected runtypes for filtering.

    Raises:
        SystemExit: If invalid metric/runtype combinations are used.
    """
    UNIVERSAL_METRICS = ["Total_Runtime_hours", "Total_Yield_Gb"]
    RUNTYPE_SPECIFIC_METRICS = ["Q30_%", "PhiX_Error_Rate_%"]

    # Calculate how many runtypes are being plotted
    # If no runtype filter specified, we're plotting all 3 types
    num_runtypes = 3 if not runtypes else len(runtypes)
    
    # Check if user requested any runtype-specific metrics
    requested_specific_metrics = [m for m in metrics if m in RUNTYPE_SPECIFIC_METRICS]

    # ERROR: User wants to plot Q30 or PhiX Error for multiple runtypes
    if num_runtypes >= 2 and requested_specific_metrics:
        logger.error("=" * 60)
        logger.error("INVALID METRIC SELECTION")
        logger.error("=" * 60)
        if num_runtypes == 3:
            logger.error("When plotting ALL run types (no --runtype filter),")
        else:
            logger.error(f"When plotting {num_runtypes} run types ({', '.join(runtypes)}),")
        logger.error("you can ONLY use these metrics:")
        logger.error(f"  - {', '.join(UNIVERSAL_METRICS)}")
        logger.error("")
        logger.error("You requested these incompatible metrics:")
        logger.error(f"  - {', '.join(requested_specific_metrics)}")
        logger.error("")
        logger.error("To use Q30_% or PhiX_Error_Rate_%, specify exactly ONE run type:")
        logger.error("  Example: --runtype RACP --metrics Q30_% PhiX_Error_Rate_%")
        logger.error("=" * 60)
        sys.exit(1)


# ============================================================================
# Plotting Helper
# ============================================================================
def plot_metric(df, metric, label_prefix=""):
    """
    Plot a scatter plot for a given metric from the DataFrame.
    
    This function handles:
    - Data type conversion (ensures numeric values)
    - NaN filtering (removes invalid/missing values)
    - Creates scatter points for the metric values across different run folders
    
    Parameters:
        df (pd.DataFrame): DataFrame with metric and RunFolder columns.
        metric (str): Column name to plot on y-axis.
        label_prefix (str): Prefix for legend label (e.g., "RACP - " when grouping by runtype).
    """
    # Convert metric column to numeric, coercing errors to NaN
    df[metric] = pd.to_numeric(df[metric], errors="coerce")
    
    # Track how many NaN values we're dropping
    initial = len(df)
    df = df[df[metric].notna()]
    dropped = initial - len(df)

    if dropped > 0:
        logger.warning(f"Dropped {dropped} NaN values for {label_prefix}{metric}")
    
    # If no valid data remains, warn and return
    if df.empty:
        logger.warning(f"No valid data to plot for {label_prefix}{metric}")
        return

    # Create scatter plot
    plt.scatter(df["RunFolder"], df[metric], alpha=0.8, label=f"{label_prefix}{metric}")
    logger.info(f"Plotted {len(df)} points for {label_prefix}{metric}")


# ============================================================================
# Add Q30 Threshold Line
# ============================================================================
def add_q30_threshold_line(runtypes, metrics):
    """
    Add a horizontal line representing the Q30 threshold if conditions are met.
    
    Different run types have different Q30 quality thresholds:
    - Panel: 80%
    - RACP: 75% 
    - Genotyping: 75% 
    
    The threshold line is only added when:
    1. Exactly one run type is being plotted (avoids confusion with mixed thresholds)
    2. Q30_% is in the selected metrics (threshold only makes sense for Q30)

    Parameters:
        runtypes (list): List of selected runtypes.
        metrics (list): List of selected metrics.
    """
    Q30_THRESHOLDS = {
        "Panel": 80,
        "RACP": 75,
        "Genotyping": 75
    }
    
    # Only add threshold if exactly one runtype and Q30 is being plotted
    if runtypes and len(runtypes) == 1 and "Q30_%" in metrics:
        runtype = runtypes[0]
        threshold = Q30_THRESHOLDS.get(runtype)
        if threshold:
            plt.axhline(y=threshold, color='red', linestyle='--', linewidth=2,
                        label=f'{runtype} Q30 Threshold ({threshold}%)', zorder=1)
            logger.info(f"Added Q30 threshold line at {threshold}% for {runtype}")


# ============================================================================
# Add Month Separators to Plot
# ============================================================================
def add_month_separators(df, ax):
    """
    Add vertical lines and labels to separate different months in the plot.
    
    When plotting data from multiple months, this creates clear visual boundaries
    between each month with:
    - Red dashed vertical lines at month boundaries
    - Month labels at the top of the plot
    
    Parameters:
        df (pd.DataFrame): DataFrame containing archive_month column
        ax (matplotlib.axes.Axes): Plot axes to add separators to
    """
    # Only add separators if we have archive_month data
    if 'archive_month' not in df.columns or df['archive_month'].isna().all():
        return
    
    # Get list of unique months (should already be sorted from data loading)
    unique_months = sorted(df['archive_month'].dropna().unique())
    
    # Don't add separators if only one month
    if len(unique_months) <= 1:
        return
    
    logger.info(f"Adding month separators for: {unique_months}")
    
    # Sort dataframe by RunFolder and get unique run folders
    df_sorted = df.sort_values("RunFolder")
    run_folders = df_sorted["RunFolder"].unique()
    
    # Find boundaries between months by checking when archive_month changes
    month_boundaries = []
    for i in range(len(run_folders) - 1):
        curr_month = df_sorted[df_sorted["RunFolder"] == run_folders[i]]["archive_month"].iloc[0]
        next_month = df_sorted[df_sorted["RunFolder"] == run_folders[i+1]]["archive_month"].iloc[0]
        if curr_month != next_month:
            # Place line between the two runs (at position i + 0.5)
            month_boundaries.append(i + 0.5)
    
    # Draw vertical lines at month boundaries
    y_min, y_max = ax.get_ylim()
    for boundary in month_boundaries:
        ax.axvline(boundary, color='gray', linestyle='--', linewidth=2, alpha=0.7)
    
    # Add month labels at the top of the plot
    month_positions = []
    for month in unique_months:
        # Find all runs from this month
        month_runs = df_sorted[df_sorted["archive_month"] == month]["RunFolder"].unique()
        if len(month_runs) > 0:
            # Find start and end positions of this month's runs
            start_idx = list(run_folders).index(month_runs[0])
            end_idx = list(run_folders).index(month_runs[-1])
            # Place label in the middle of the month's runs
            mid_pos = (start_idx + end_idx) / 2
            month_positions.append((mid_pos, month))
    
    # Draw text labels for each month
    for pos, month in month_positions:
        ax.text(pos, y_max * 1.05, f"Month: {month}", 
               ha='center', va='bottom', fontsize=10, fontweight='bold',
               bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.7))


# ============================================================================
# Main Function
# ============================================================================
def main():
    """
    Main orchestration function that:
    1. Parses command line arguments
    2. Handles --list-months flag
    3. Loads data from archived months
    4. Validates metric/runtype compatibility
    5. Creates and saves the plot with selected metrics
    """
    logger.info("=" * 60)
    logger.info("Starting plot_trends.py")
    logger.info("=" * 60)

    args = parse_args()
    
    # Handle --list-months flag (early exit)
    if args.list_months:
        list_available_months()
        return

    # Resolve metric aliases (e.g., "q30" → "Q30_%")
    resolved_metrics = []
    for m in args.metrics:
        key = m.lower()
        resolved_metrics.append(METRIC_ALIASES.get(key, m))
    args.metrics = resolved_metrics

    # Prepare arguments
    runtypes = args.runtype or []

    # Log inputs for debugging
    logger.info(f"Metrics requested: {', '.join(args.metrics)}")
    logger.info(f"Run type filter: {', '.join(runtypes) if runtypes else 'ALL'}")
    logger.info(f"Output file: {args.out}")
    logger.info("=" * 60)

    # Validate metric/runtype compatibility (enforces rules about which metrics work with which runtypes)
    validate_metric_compatibility(args.metrics, runtypes)

    # Parse and validate month list
    months = [m.strip() for m in args.months.split(",")]
    
    # Validate month format (must be YYYY-MM)
    import re
    for month in months:
        try:
            datetime.strptime(month, "%Y-%m")
        except ValueError:
            logger.error(f"Invalid month format: {month}")
            raise SystemExit(f"Error: Month '{month}' is not in YYYY-MM format")
    
    # Load data from archived months
    logger.info(f"Loading data from archived months: {', '.join(months)}")
    df = combine_monthly_archives(months)
    
    # Apply runtype filter if specified
    if runtypes and "RunType" in df.columns:
        df = df[df["RunType"].isin(runtypes)]
        logger.info(f"Filtered to {len(df)} rows for run types: {', '.join(runtypes)}")
    
    # Validate structure (ensure required columns exist)
    if "RunFolder" not in df.columns:
        logger.error("CSV missing required 'RunFolder' column")
        raise ValueError("CSV must contain a 'RunFolder' column for the x-axis.")
    
    # Clean RunType column if it exists
    if "RunType" in df.columns:
        df["RunType"] = df["RunType"].astype(str).str.strip()

    # Validate that requested metrics exist in CSV
    missing = [m for m in args.metrics if m not in df.columns]
    if missing:
        raise ValueError(f"Metric(s) not found: {', '.join(missing)}. "
                         f"Available columns: {', '.join(df.columns)}")

    # Set up plot_output directory inside month_0_output folder
    from pathlib import Path
    project_root = Path(__file__).resolve().parent.parent.parent.parent.parent  
    month_0_output_dir = project_root / "month_0_output"
    plot_output_dir = month_0_output_dir / "plot_output"
    plot_output_dir.mkdir(parents=True, exist_ok=True)

    # Create plot
    plt.figure(figsize=(30, 11))

    # If no runtype filter, group by RunType for separate colors/labels
    if "RunType" in df.columns and not runtypes:
        for runtype, group_df in df.groupby("RunType"):
            for metric in args.metrics:
                plot_metric(group_df.copy(), metric, label_prefix=f"{runtype} - ")
    else:
        # Filtered by runtype or no RunType column
        filtered_df = df.copy()
        for metric in args.metrics:
            plot_metric(filtered_df.copy(), metric)

    # Get current axes for adding separators
    ax = plt.gca()
    
    # Add month separators if multiple months present
    add_month_separators(df, ax)
    
    # Add Q30 threshold line if appropriate
    add_q30_threshold_line(runtypes, args.metrics)

    # Final plot settings
    plt.grid(True, linestyle='--', alpha=0.4)
    
    # Build title with month info if applicable
    title = f"{args.metrics} by Run"
    if runtypes:
        title += f" - {', '.join(runtypes)}"
    if len(months) > 1:
        title += f" ({len(months)} months)"
    
    plt.title(title)
    plt.xlabel("RunFolder")
    plt.ylabel("Metric Value")
    plt.xticks(rotation=60, ha="right")
    
    # Show legend if we have multiple series
    if "RunType" in df.columns or len(months) > 1:
        plt.legend()
    
    plt.ylim(bottom=0)  # Start y-axis at 0
    plt.tight_layout()
    
    # Save plot
    plot_path = plot_output_dir / args.out
    plt.savefig(plot_path, dpi=150)
    logger.info(f"Plot successfully saved: {plot_path}")

    logger.info("=" * 60)
    logger.info("plot_trends.py completed successfully")
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
    