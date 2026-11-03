"""
Horizontal Coverage Analysis - Exon Coverage QC Reporting

Analyses exon coverage data from Panel sequencing runs to identify samples with 
insufficient coverage. Generates boxplots showing coverage ≥20x distribution 
per run and exports outlier samples falling below the 90% threshold.
"""

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

# Set up logging
from utils.logger import create_logger
from utils.config import OUTPUT_ROOT, RESULTS_DIR
logger = create_logger('HorCov.log')

def parse_arguments():
    """
    Parse command-line arguments for month and year selection.
    """
    parser = argparse.ArgumentParser(
        description="Analyse exon coverage data from sequencing runs."
    )
    parser.add_argument(
        "--months",
        type=str,
        help="Comma-separated list of YYYY-MM dates to analyze (e.g., '2025-01,2025-05,2026-01')"
    )
    return parser.parse_args()


def parse_month_specifications(month_spec):
    """
    Parse month specification string into list of (year, month) tuples.
    
    Args:
        month_spec: String like "2025-01,2025-05,2026-01"
    
    Returns:
        List of (year, month) tuples
    """
    if not month_spec:
        return None
    
    months = []
    for spec in month_spec.split(','):
        spec = spec.strip()
        try:
            year, month = spec.split('-')
            year = int(year)
            month = int(month)
            if month < 1 or month > 12:
                logger.error(f"Invalid month in '{spec}'. Month must be 1-12.")
                sys.exit(f"Error: Invalid month in '{spec}'. Month must be 1-12.")
            months.append((year, month))
        except ValueError:
            logger.error(f"Invalid format '{spec}'. Use YYYY-MM format.")
            sys.exit(f"Error: Invalid format '{spec}'. Use YYYY-MM format.")
    
    return months


def extract_run_date_from_interop(btch_dir):
    """
    Extract run date from interop file in BTCH directory.
    Looks for files like: 250101_MN00000_0001_A000ABCDEF_interop...
    
    Args:
        btch_dir: Path to BTCH directory
    
    Returns:
        (year, month, day) tuple or None if not found
    """
    # Look for interop files in the BTCH directory
    interop_files = list(btch_dir.glob("*_interop*"))
    
    for interop_file in interop_files:
        filename = interop_file.name
        # Extract date from filename (YYMMDD at start)
        date_match = re.match(r'^(\d{6})_', filename)
        if date_match:
            date_str = date_match.group(1)
            try:
                # Parse YYMMDD format
                year = int("20" + date_str[0:2])  # Convert YY to YYYY
                month = int(date_str[2:4])
                day = int(date_str[4:6])
                return (year, month, day)
            except ValueError:
                continue
    
    return None


def get_btch_directories(base_years):
    """
    Find all BTCH directories across specified years.
    
    Args:
        base_years: Set of years to search
    
    Returns:
        List of Path objects for BTCH directories
    """
    base_dir = RESULTS_DIR
    btch_dirs = []
    
    for year in base_years:
        year_dir = base_dir / str(year)
        if year_dir.exists():
            # Find all directories containing 'BTCH' in their name
            btch_dirs.extend(year_dir.glob("*BTCH*"))
    
    return btch_dirs


def filter_btch_by_months(btch_dirs, target_months):
    """
    Filter BTCH directories by extracting run dates from interop files.
    
    Args:
        btch_dirs: List of BTCH directory paths
        target_months: List of (year, month) tuples to include, or None for all
    
    Returns:
        Dictionary mapping btch_dir to (year, month, day) tuple
    """
    btch_with_dates = {}
    
    for btch_dir in btch_dirs:
        run_date = extract_run_date_from_interop(btch_dir)
        if run_date:
            year, month, day = run_date
            # Check if this matches target months (if specified)
            if target_months is None or (year, month) in set(target_months):
                btch_with_dates[btch_dir] = run_date
                logger.info(f"Including {btch_dir.name} (run date: {year}-{month:02d}-{day:02d})")
    
    return btch_with_dates


def extract_btch_id(path_str):
    """
    Extract BTCH ID from path, handling both formats:
    - 250519_BTCH600744_1
    - BTCH600772_2
    
    Args:
        path_str: String representation of file path
    
    Returns:
        BTCH ID string or "unknown"
    """
    # Try to match BTCH pattern with optional date prefix
    btch_match = re.search(r'(?:\d{6}_)?(BTCH\d+(?:_\d+)?)', path_str)
    if btch_match:
        return btch_match.group(1)
    return "unknown"


def main():
    """Main execution function."""
    logger.info("=== Starting Horizontal Coverage Analysis ===")
    
    # Parse command-line arguments
    args = parse_arguments()
    
    # Parse month specifications
    target_months = parse_month_specifications(args.months)
    
    # Determine which years to search based on target months
    if target_months:
        base_years = set(year for year, month in target_months)
        logger.info(f"Searching for data in years: {sorted(base_years)}")
        logger.info(f"Target months: {', '.join(f'{y}-{m:02d}' for y, m in target_months)}")
    else:
        # If no months specified, use current year
        base_years = {datetime.now().year}
        logger.info(f"No months specified. Searching all data in year: {datetime.now().year}")
    
    # Find all BTCH directories
    btch_dirs = get_btch_directories(base_years)
    logger.info(f"Found {len(btch_dirs)} BTCH directories")
    
    # Filter BTCH directories by month and get dates
    btch_with_dates = filter_btch_by_months(btch_dirs, target_months)
    logger.info(f"After filtering by month: {len(btch_with_dates)} BTCH directories")
    
    if not btch_with_dates:
        logger.error("No BTCH directories found matching the specified criteria.")
        sys.exit("No BTCH directories found matching the specified criteria.")
    
    # Create timestamp for output files
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Create month suffix for filenames if months specified
    month_suffix = ""
    if target_months:
        month_strs = [f"{y}{m:02d}" for y, m in sorted(target_months)]
        month_suffix = f"_months_{'_'.join(month_strs)}"
    
    # Defines the output directory
    output_dir = OUTPUT_ROOT / "month_1_output_HorCov"
    output_dir.mkdir(parents=True, exist_ok=True)  # Create output directory if it doesn't exist
    
    # --- 1) Find all W*_exoncoverage_bam.txt files in filtered BTCH directories ---
    coverage_files = []
    for btch_dir in btch_with_dates.keys():
        coverage_files.extend(btch_dir.glob("W*/W*_exoncoverage_bam.txt"))  # Find all files matching the pattern
    
    if not coverage_files:
        logger.error("No W*_exoncoverage_bam.txt files found under the expected pattern.")
        sys.exit("No W*_exoncoverage_bam.txt files found under the expected pattern.")  # Exit if no files found under the expected pattern
    
    logger.info(f"Found {len(coverage_files)} coverage files to process")
    
    coverage_data = []  # Define an empty list to store results
    
    # Iterate over each coverage report file
    for report_path in coverage_files:
        sample_id = report_path.parent.name  # Extract the W number (sample ID) from the folder name
        run_id = extract_btch_id(str(report_path))  # Extract the BTCH ID from the file path
        
        # --- 2) Parse the coverage report and compute coverage weighted %≥20x ---
        try:
            coverage_df = pd.read_csv(
                report_path, sep=r"\s+", engine="python", comment="#",  # Use regex for whitespace separation, Python engine, and ignore comment lines
                names=["EXON", "LENGTH", "COVERED", "PERCENT"], usecols=[0, 1, 2, 3]  # Specify column names and use only 4 relevant columns
            ).dropna()  # Drop rows with any NaN values
        except Exception as e:
            logger.warning(f"Could not read {report_path}: {e}")
            continue
        
        # Sum exon lengths and covered bases at ≥20x
        total_length = coverage_df["LENGTH"].astype(int).sum()  # Total exon length
        total_covered = coverage_df["COVERED"].astype(int).sum()  # Total bases covered at ≥20x
        pct_20x = (100.0 * total_covered / total_length)  # Calculate percentage of bases covered at ≥20x
        
        # Append results to the list
        coverage_data.append({
            "Run_ID": run_id,
            "Sample": sample_id,
            "pct_20x": pct_20x,
            "report_path": str(report_path)
        })
    
    # Create a DataFrame from all parsed data and save to CSV
    all_samples_df = pd.DataFrame(coverage_data).sort_values(["Run_ID", "Sample"])
    coverage_csv = output_dir / f"coverage_all_samples{month_suffix}_{timestamp}.csv"
    all_samples_df.to_csv(coverage_csv, index=False)
    logger.info(f"Saved coverage data → {coverage_csv}")
    
    # --- 3) Identify and save outliers elsewhere (here defined as <90% at ≥20x) ---
    coverage_threshold = 90.0  # Define the coverage threshold for outliers
    outliers_samples_df = all_samples_df[all_samples_df["pct_20x"] < coverage_threshold].copy()  # Filter out true outliers
    outliers_csv = output_dir / f"outliers_below_{int(coverage_threshold)}pct{month_suffix}_{timestamp}.csv"
    outliers_samples_df.to_csv(outliers_csv, index=False)  # Save outliers to separate CSV
    logger.info(f"Identified {len(outliers_samples_df)} outliers (<{coverage_threshold}% coverage) → {outliers_csv}")
    
    # --- 4) Plot per-run boxplots with overlay dots (≥90%) ---
    # Create mapping of Run_ID to run date for x-axis labels
    run_id_to_date = {}
    for btch_dir, run_date in btch_with_dates.items():
        run_id = extract_btch_id(str(btch_dir))
        if run_id != "unknown":
            year, month, day = run_date
            run_id_to_date[run_id] = f"{year}-{month:02d}-{day:02d}"
    
    run_ids = all_samples_df["Run_ID"].dropna().unique().tolist()  # Get unique Run_ID values
    if not run_ids:
        logger.error("No Run_ID values found")
        sys.exit("No Run_ID values found.")  # Exit if no Run_ID values found
    
    coverage_per_run = []  # List to hold coverage data per run
    for run_id in run_ids:
        run_values = all_samples_df.loc[all_samples_df["Run_ID"] == run_id, "pct_20x"]  # Get pct_20x values for the current run
        coverage_per_run.append(run_values[run_values >= coverage_threshold].to_numpy())  # Append values ≥90% to the list
    
    # Figure size calculation with more generous width per run
    fig_width = max(12, 0.6 * len(run_ids))  # At least 12 inches wide, or 0.6 inches per run
    fig, ax = plt.subplots(figsize=(fig_width, 10))
    
    # Spacing between boxplots
    x_positions = np.arange(len(run_ids)) * 0.5  # Increased spacing from 0.4 to 0.5
    
    boxplot_result = ax.boxplot(
        coverage_per_run,
        vert=True,  # Vertical boxplots
        positions=x_positions,  # Set spacing
        widths=0.3,  # Slightly wider boxes
        patch_artist=True,  # Enable custom box colors
        showmeans=True,  # Display mean
        meanline=True,  # Use a line for the mean
        labels=run_ids,  # Run names on x-axis
        meanprops=dict(color="black", linewidth=0.8, linestyle="dashed"),  # Mean line properties
        medianprops=dict(color="black", linewidth=1.4),  # Median line properties
        boxprops=dict(linewidth=1.0),
        whiskerprops=dict(linewidth=1.0),
        capprops=dict(linewidth=1.0),
        showfliers=False,  # Hide statistical outliers
    )
    
    # Apply colours to boxplots
    cmap = cm.get_cmap("tab20", len(run_ids))  # Colormap with distinct colors
    for i, patch in enumerate(boxplot_result["boxes"]):
        patch.set_facecolor(cmap(i))  # Set box color
        patch.set_alpha(0.5)  # Set transparency
    
    # Overlay sample points
    for x_pos, run_id in zip(x_positions, run_ids):
        run_sample_dfs = all_samples_df[
            (all_samples_df["Run_ID"] == run_id) &
            (all_samples_df["pct_20x"] >= coverage_threshold)]
        
        for _, sample in run_sample_dfs.iterrows():
            ax.scatter(x_pos, sample["pct_20x"], s=14, alpha=0.6, edgecolors="none", c="grey")
    
    # Axes settings
    ax.set_ylim(95, 100)
    ax.set_ylabel("Coverage ≥20× (% of target bases per sample)", fontsize=11)
    ax.set_xlabel("Sequencing Run", fontsize=11)
    
    # Update x-axis labels to show date alongside BTCH ID 
    x_labels = []
    for run_id in run_ids:
        if run_id in run_id_to_date:
            x_labels.append(f"{run_id} ({run_id_to_date[run_id]})")
        else:
            x_labels.append(run_id)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, fontsize=8)  # Smaller font for better fit
    
    # Short, clear title
    title_suffix = ""
    if target_months:
        month_strs = [f"{y}-{m:02d}" for y, m in sorted(target_months)]
        if len(month_strs) <= 3:
            title_suffix = f" ({', '.join(month_strs)})"
        else:
            title_suffix = f" ({len(month_strs)} months)"
    
    ax.set_title(f"Coverage ≥20x for Samples with ≥90% Target Coverage{title_suffix}", fontsize=12, pad=15)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")  # 45 degree angle instead of 60
    ax.grid(True, axis="y", linestyle=":", alpha=0.8)
    ax.grid(True, axis="x", linestyle=":", alpha=0.8)
    
    # Add legend for mean and median lines
    legend_elements = [
        Line2D([0], [0], color='black', lw=1.4, label='Median (solid)'),
        Line2D([0], [0], color='black', lw=0.8, linestyle='--', label='Mean (dashed)')
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9)
    
    # Tight layout to prevent cutoff
    fig.tight_layout(pad=2.0)  # Extra padding to prevent label cutoff
    
    # Clear filename 
    out_png = output_dir / f"coverage_20x{month_suffix}_{timestamp}.png"
    fig.savefig(out_png, dpi=160, bbox_inches='tight')  # bbox_inches='tight' ensures nothing gets cut off
    plt.close(fig)
    logger.info(f"Boxplot saved → {out_png}")
    logger.info("=== Analysis complete ===")


if __name__ == "__main__":
    main()
    