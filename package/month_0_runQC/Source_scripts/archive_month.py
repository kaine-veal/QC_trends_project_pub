"""
archive_month.py

This script processes raw Illumina sequencing run data for a specific month and creates
a monthly archive CSV. This is the first step in the QC workflow - archive monthly data,
then use plot_trends.py to visualise it.

Key features
1. Processes only runs from the specified month
2. Parses run dates from folder names
3. Determines run types (RACP, Genotyping, Panel) from SampleSheet.csv
4. Runs Illumina's interop summary tool to extract metrics
5. Calculates total runtime by parsing timestamps
6. Saves results to monthly archive directory (See readme for structure)
"""

import subprocess  # To run external shell commands from Python (e.g. interop summary)
import csv  # To write results into a CSV file
import re  # To search for patterns in strings (e.g. timestamps)
import argparse # To handle command line arguments
from pathlib import Path  # To handle file and directory paths
from datetime import datetime, timezone  # This is for extracting the TAT for each run

# Set up logging
from utils.logger import create_logger
from utils.config import OUTPUT_ROOT
logger = create_logger('run_QC.log')

# Archive location for monthly CSVs
ARCHIVE_BASE = OUTPUT_ROOT / "qc_archive_runqc"


# ============================================================================
# Run Date Parsing
# ============================================================================
def get_run_date(name: str):
    """
    The purpose of this function is to get the run date from the run file name to enable
    correct filtering of runs in run directory by month.
    
    Supports two formats:
    - YYMMDD (e.g., 250101_MN00000_0001_A000ABCDEF)
    - YYYYMMDD (e.g., 20251117_...)
    
    Parameters:
        name (str): Run folder name (e.g., '250101_MN00000_0001_A000ABCDEF')
    
    Returns:
        datetime: Parsed date object, or None if parsing fails
    """
    s = str(name)

    # Match YYMMDD format first (6 digits at start)
    m = re.match(r"^(\d{2})(\d{2})(\d{2})", s)
    if m:
        try:
            yy, mm, dd = m.groups()
            year = 2000 + int(yy)  
            parsed_date = datetime(year, int(mm), int(dd))
            logger.debug(f"Parsed date from '{name}' as YYMMDD: {parsed_date.strftime('%Y-%m-%d')}")
            return parsed_date
        except ValueError:
            pass

    # Match YYYYMMDD format (8 digits)
    m = re.search(r"\b(\d{8})\b", s)
    if m:
        try:
            parsed_date = datetime.strptime(m.group(1), "%Y%m%d")
            logger.debug(f"Parsed date from '{name}' as YYYYMMDD: {parsed_date.strftime('%Y-%m-%d')}")
            return parsed_date
        except ValueError:
            pass

    logger.warning(f"Could not parse date from run folder name: {name}")
    return None


# ============================================================================
# Run Type Determination
# ============================================================================
def determine_run_type(run_folder: Path) -> str:
    """
    Determine the run type by checking the SampleSheet.csv file for keywords.
    
    The function searches for specific keywords in the sample sheet:
    - "RACP" → RACP sequencing run
    - "Genotyping" → Genotyping run
    - "Panel" → Panel sequencing run
    
    The search is case-insensitive and checks multiple variations of each keyword.
    Priority order: RACP > Panel > Genotyping (because some sheets may contain multiple keywords)
    
    Parameters:
        run_folder (Path): Path to the run folder
    
    Returns:
        str: One of 'RACP', 'Genotyping', 'Panel', or 'Unknown'
    """
    runtype = "Unknown"  # Default value in case the sample sheet is missing or doesn't contain run name
    sample_sheet_path = run_folder / "SampleSheet.csv"  # Define the expected path to the sample sheet for this run
    
    if sample_sheet_path.exists():
        try:
            with open(sample_sheet_path, "r", encoding="utf-8") as ss:  # Open the file in read mode with UTF-8 encoding
                content = ss.read().upper()  # Read entire content and convert to uppercase for case-insensitive matching
                
                # Check in priority order (RACP > Panel > Genotyping)
                # Use case-insensitive matching by converting content to uppercase
                if "RACP" in content:
                    runtype = "RACP"
                elif "PANEL" in content:
                    runtype = "Panel"
                elif "GENOTYPING" in content or "GENOTYPE" in content:
                    runtype = "Genotyping"
                
            logger.debug(f"Determined run type for {run_folder.name}: {runtype}")
        except Exception as e:
            logger.error(f"Could not read SampleSheet for {run_folder.name}: {e}")
    else:
        logger.debug(f"SampleSheet.csv not found for {run_folder.name}")
    
    return runtype


# ============================================================================
# Runtime Calculation
# ============================================================================
def get_runtime_hours(run_folder: Path) -> float | str:
    """
    The purpose of this function is to calculate the total runtime of RACP, Genotyping, and Panel runs in hours.
    It does this by finding the start time (earliest timestamp from any RunSetup.imf1 file in the 
    InstrumentAnalyticsLogs directory) and the end time (timestamp in the RTA Complete.txt file) of each run, 
    then calculating the difference between them in hours.
    
    Parameters:
        run_folder (Path): Path to the run folder
    
    Returns:
        float: Runtime in hours (rounded to 2 decimal places)
        str: "NA" if runtime cannot be calculated
    """
    logger.debug(f"Calculating runtime for {run_folder.name}")
    
    # Start: timestamp in InstrumentAnalyticsLogs/*RunSetup.imf1
    start_time = None  # Default to None if no RunSetup.imf1 files are found
    instrumentanalyticslogs = run_folder / "InstrumentAnalyticsLogs"  # Path to Instrument Analytics Logs folder

    if instrumentanalyticslogs.is_dir():  # Only proceed if the InstrumentAnalyticsLogs directory exists
        timestamps = []  # List to store all parsed datetime objects found

        for path in instrumentanalyticslogs.rglob("*RunSetup.imf1"):  # Recursively search for files ending with 'RunSetup.imf1'
            timestamp_match = re.search(r'(\d{4}-\d{2}-\d{2}-\d{2}\.\d{2}\.\d{2})', path.name)  # Extracts timestamp from filename using regex (pattern-matching language), matching the format YYYY-MM-DD-HH.MM.SS

            if timestamp_match:
                timestamps.append(datetime.strptime(timestamp_match.group(1), "%Y-%m-%d-%H.%M.%S").replace(tzinfo=timezone.utc))  # Converts matched timestamp string into a datetime object and assigns UTC timezone
                
        if timestamps: 
            start_time = min(timestamps)  # If any timestamps were found, set start_time to the earliest one
            logger.debug(f"Found start time: {start_time}")
        else:
            logger.warning(f"No RunSetup.imf1 files found for {run_folder.name}")
    else:
        logger.warning(f"InstrumentAnalyticsLogs directory not found for {run_folder.name}")

    # End: timestamp in RTAComplete.txt
    end_time = None  # Default to None if no RTAComplete.txt file is found
    rta_complete = run_folder / "RTAComplete.txt"  # Path to RTA Complete file

    if rta_complete.exists():
        rta_contents = rta_complete.read_text(errors="ignore").strip()  # ensures rta_contents is a string with no leading/trailing whitespace ready for regex searching

        # Format A: '6/19/2025,03:19:21.929,Illumina RTA 1.18.54.4'
        match_a = re.search(r'(\d{1,2})/(\d{1,2})/(\d{4}),(\d{2}):(\d{2}):(\d{2})', rta_contents)
        if match_a:
            MM, DD, YYYY, hh, mm, ss = match_a.groups()
            end_time = datetime(int(YYYY), int(MM), int(DD), int(hh), int(mm), int(ss), tzinfo=timezone.utc)
            logger.debug(f"Found end time (Format A): {end_time}")

        else:
            # Format B: 'RTA 2.11.4.0 completed on 3/11/2025 2:04:22 PM'
            match_b = re.search(r'(\d{1,2})/(\d{1,2})/(\d{4})\s+(\d{1,2}):(\d{2}):(\d{2})\s*(AM|PM)', rta_contents, re.I)
            if match_b:
                MM, DD, YYYY, hh, mm, ss, ampm = match_b.groups()
                hh = int(hh) % 12 + (12 if ampm.upper() == "PM" else 0)
                end_time = datetime(int(YYYY), int(MM), int(DD), hh, int(mm), int(ss), tzinfo=timezone.utc)
                logger.debug(f"Found end time (Format B): {end_time}")
            else:
                logger.warning(f"Could not parse RTAComplete.txt timestamp for {run_folder.name}")
    else:
        logger.warning(f"RTAComplete.txt not found for {run_folder.name}")

    if start_time and end_time:  # Calculate runtime only if both timestamps are available
        runtime = round((end_time - start_time).total_seconds() / 3600, 2)  # Convert total seconds to hours and round to 2 decimal places
        logger.info(f"Runtime for {run_folder.name}: {runtime} hours")
        return runtime
    
    logger.warning(f"Could not calculate runtime for {run_folder.name} - missing start or end time")
    return "NA"


# ============================================================================
# InterOp Metrics Extraction
# ============================================================================
def extract_interop_metrics(run_folder: Path) -> dict | None:
    """
    Run Illumina interop summary and extract key quality metrics.
    
    This function:
    1. Executes the 'interop summary' command on the run folder
    2. Parses the output to find the "Total" line with aggregate metrics
    3. Extracts Q30%, PhiX error rate, and total yield
    
    Parameters:
        run_folder (Path): Path to the run folder
    
    Returns:
        dict: Metrics dictionary with keys:
            - Q30_%: Percentage of reads with Q30 or higher
            - PhiX_Error_Rate_%: PhiX error rate percentage
            - Total_Yield_Gb: Total yield in gigabases
        None: If extraction fails
    """
    logger.debug(f"Extracting interop metrics for {run_folder.name}")
    
    try:
        # "source ~/.bashrc" - loads the server environment so 'interop' alias becomes available
        # "interop summary '{run_folder}" - runs Illumina interop summary on specific run
        cmd = f"source ~/.bashrc && interop summary '{run_folder}'"
        logger.debug(f"Executing command: {cmd}")

        # Run the command inside a Bash shell, capture both stdout and stderr
        result = subprocess.run(
            cmd,
            shell=True, executable="/bin/bash",  # Run inside Bash to use aliases
            stdout=subprocess.PIPE,  # Capture normal output to use in python
            stderr=subprocess.PIPE,  # Capture error messages
            universal_newlines=True,  # Return output as text instead of bytes
            check=True  # Raise error if command fails
        )
        logger.debug(f"interop summary completed successfully for {run_folder.name}")

    # If interop summary fails, log error and return None (run will be skipped)
    except subprocess.CalledProcessError as e:
        logger.error(f"interop summary failed for {run_folder.name}: {e}")
        logger.error(f"STDERR: {e.stderr}")
        return None

    # Takes entire output from result.stdout (one long string) and splits it into individual lines
    lines = result.stdout.splitlines()

    # Loops through each output line to find the one that starts with "Total"
    total_line = None
    for line in lines:
        if line.startswith("Total"):
            total_line = [x.strip() for x in line.split(",")]  # Split by commas and strip whitespace
            logger.debug(f"Found Total line: {total_line}")
            break

    # If total line not found, log a warning and return None
    if not total_line or len(total_line) < 7:  # Less than 7 ensures necessary metrics are present in total_line
        logger.warning(f"Could not parse metrics for {run_folder.name} - Total line missing or incomplete")
        return None

    # Extract key metrics from known column positions
    try:
        metrics = {
            "Q30_%": total_line[6],              # Column 7 = % reads ≥ Q30
            "PhiX_Error_Rate_%": total_line[4],  # Column 5 = PhiX error rate (%)
            "Total_Yield_Gb": total_line[1]      # Column 2 = Yield in Gb
        }
        logger.debug(f"Extracted metrics - Q30: {metrics['Q30_%']}%, PhiX Error: {metrics['PhiX_Error_Rate_%']}%, Yield: {metrics['Total_Yield_Gb']} Gb")
        return metrics
    except IndexError:
        logger.error(f"Parsing issue for {run_folder.name} - column index out of range")
        return None


# ============================================================================
# Month Validation
# ============================================================================
def validate_month_format(month: str) -> bool:
    """
    Validate that month string is in YYYY-MM format.
    
    Parameters:
        month (str): Month string to validate
    
    Returns:
        bool: True if valid format, False otherwise
    """
    try:
        datetime.strptime(month, "%Y-%m")
        return True
    except ValueError:
        return False


def check_month_already_archived(month: str) -> bool:
    """
    Check if a month has already been archived.
    
    Parameters:
        month (str): Month in YYYY-MM format
    
    Returns:
        bool: True if archive exists, False otherwise
    """
    month_dir = ARCHIVE_BASE / month
    csv_path = month_dir / f"interop_metrics_{month}.csv"
    return csv_path.exists()


# ============================================================================
# Main Processing Function
# ============================================================================
def process_month(run_directory: Path, month: str, runtype_filter: list = None) -> None:
    """
    Process all runs from a specific month and create archive CSV.
    
    This is the main orchestration function that:
    1. Filters runs to those in the specified month
    2. Processes each run to extract metrics
    3. Writes results to monthly archive CSV
    
    Parameters:
        run_directory (Path): Path to directory containing run folders
        month (str): Month to process in YYYY-MM format
        runtype_filter (list, optional): List of run types to include (e.g., ['RACP', 'Panel'])
    """
    logger.info("=" * 60)
    logger.info(f"Processing month: {month}")
    logger.info(f"Run directory: {run_directory}")
    if runtype_filter:
        logger.info(f"Run type filter: {', '.join(runtype_filter)}")
    else:
        logger.info("Run type filter: ALL")
    logger.info("=" * 60)

    # Parse target month for filtering
    target_year, target_month = map(int, month.split("-"))
    
    # Column names for first row of CSV file
    header = ["RunType", "RunFolder", "Q30_%", "PhiX_Error_Rate_%", "Total_Runtime_hours", "Total_Yield_Gb"]
    
    # Empty list to store all rows of results before writing them to CSV
    rows = []
    
    # Count total folders for progress tracking
    all_folders = [f for f in run_directory.iterdir() if f.is_dir()]
    logger.info(f"Found {len(all_folders)} total folders in run directory")
    
    # Counters for reporting
    processed_count = 0
    skipped_count = 0
    wrong_month_count = 0
    runtype_filtered_count = 0
    unknown_runtype_count = 0  # Track runs with Unknown runtype
    
    # Filter to directories starting with 6-digit date format (YYMMDD) only
    run_folders = [f for f in run_directory.iterdir() 
                   if f.is_dir() and re.match(r"^\d{6}_", f.name)]
    
    # Iterate through each folder in the run directory and sort chronologically
    for run_folder in sorted(run_folders, key=lambda p: get_run_date(p.name) or datetime.min):
        
        # Parse the run date from folder name
        rd = get_run_date(run_folder.name)
        
        # Skip runs not in target month
        if not rd or rd.year != target_year or rd.month != target_month:
            logger.debug(f"Skipping {run_folder.name} - not in target month {month}")
            wrong_month_count += 1
            continue
        
        # Skip runs without InterOp folder
        if not (run_folder / "InterOp").exists():
            logger.warning(f"Skipping {run_folder.name} - no InterOp folder found")
            skipped_count += 1
            continue
        
        # Determine run type from SampleSheet.csv
        runtype = determine_run_type(run_folder)
        
        # Skip runs with Unknown run type (no valid keyword found in SampleSheet)
        if runtype == "Unknown":
            logger.warning(f"Skipping {run_folder.name} - run type 'Unknown' (no RACP/Genotyping/Panel keyword found in SampleSheet)")
            unknown_runtype_count += 1
            continue
        
        # Apply runtype filter if specified (skip runs not in filter)
        if runtype_filter and runtype not in runtype_filter:
            logger.debug(f"Skipping {run_folder.name} - run type '{runtype}' not in filter {runtype_filter}")
            runtype_filtered_count += 1
            continue
        
        # Show which run is currently being processed
        logger.info(f"Processing run: {run_folder.name} (RunType: {runtype})")
        
        # Extract metrics using interop summary
        metrics = extract_interop_metrics(run_folder)
        if metrics is None:
            skipped_count += 1
            continue
        
        # Calculate runtime in hours
        runtime_hours = get_runtime_hours(run_folder)
        
        # Store extracted metrics in rows list for CSV writing
        rows.append([
            runtype,
            run_folder.name,
            metrics["Q30_%"],
            metrics["PhiX_Error_Rate_%"],
            runtime_hours,
            metrics["Total_Yield_Gb"]
        ])
        logger.info(f"Successfully processed {run_folder.name} - added to results")
        processed_count += 1
    
    # Create archive directory if it doesn't exist
    month_dir = ARCHIVE_BASE / month
    month_dir.mkdir(parents=True, exist_ok=True)
    
    # Define output CSV path
    output_csv = month_dir / f"interop_metrics_{month}.csv"
    
    # Write results to CSV
    logger.info(f"Writing {len(rows)} rows to CSV: {output_csv}")
    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)  # Write column headers as first row
        writer.writerows(rows)   # Write one row per sequencing run
    
    logger.info(f"Metrics successfully saved to {output_csv}")
    
    # Print summary
    logger.info("=" * 60)
    logger.info("Processing Summary:")
    logger.info(f"  - Total folders found: {len(all_folders)}")
    logger.info(f"  - Wrong month: {wrong_month_count}")
    logger.info(f"  - Successfully processed: {processed_count}")
    logger.info(f"  - Unknown runtype (skipped): {unknown_runtype_count}")
    logger.info(f"  - Filtered by runtype: {runtype_filtered_count}")
    logger.info(f"  - Skipped (errors/missing data): {skipped_count}")
    logger.info("=" * 60)
    logger.info(f"Archive created: {output_csv}")
    logger.info("=" * 60)


# ============================================================================
# Command Line Interface
# ============================================================================
def parse_args():
    """
    Parse command line arguments for the script.
    
    Returns:
        argparse.Namespace: Parsed arguments containing:
            - rundir: Path to directory with sequencing run folders
            - month: Month to archive (YYYY-MM format)
            - runtype: Optional list of run types to include
            - force: Whether to overwrite existing archive
    """
    ap = argparse.ArgumentParser(
        description="Archive Illumina sequencing run metrics for a specific month",
    )
    ap.add_argument("--rundir", required=True, help="Directory containing run folders")
    ap.add_argument("--month", required=True, help="Month to archive in YYYY-MM format (e.g., 2025-01)")
    ap.add_argument("--runtype", nargs="*", choices=["RACP", "Genotyping", "Panel"],
                    help="Filter by run type(s). If not specified, all run types are included.")
    ap.add_argument("--force", action="store_true", 
                    help="Overwrite existing archive if it exists")
    
    return ap.parse_args()


# ============================================================================
# Main Entry Point
# ============================================================================
def main():
    """
    Main entry point for the script.
    
    Handles argument parsing, validation, and orchestrates the archiving process.
    """
    args = parse_args()
    
    logger.info("=" * 60)
    logger.info("Starting archive_month.py")
    logger.info("=" * 60)
    
    # Validate month format
    if not validate_month_format(args.month):
        logger.error(f"Invalid month format: {args.month}")
        logger.error("Month must be in YYYY-MM format (e.g., 2025-01)")
        raise SystemExit("Error: Invalid month format")
    
    # Check if month is in the future
    month_date = datetime.strptime(args.month, "%Y-%m")
    if month_date > datetime.now():
        logger.error(f"Month {args.month} is in the future")
        raise SystemExit("Error: Cannot archive future months")
    
    # Check if already archived
    if check_month_already_archived(args.month) and not args.force:
        logger.error(f"Month {args.month} has already been archived")
        logger.error("Use --force to overwrite existing archive")
        csv_path = ARCHIVE_BASE / args.month / f"interop_metrics_{args.month}.csv"
        logger.error(f"Existing archive: {csv_path}")
        raise SystemExit("Error: Month already archived (use --force to overwrite)")
    
    # Validate run directory
    run_directory = Path(args.rundir)
    if not run_directory.exists():
        logger.error(f"Run directory does not exist: {run_directory}")
        raise FileNotFoundError(f"Run directory not found: {run_directory}")
    
    if not run_directory.is_dir():
        logger.error(f"Run directory path is not a directory: {run_directory}")
        raise NotADirectoryError(f"Not a directory: {run_directory}")
    
    # Process the month
    try:
        process_month(run_directory, args.month, args.runtype)
        logger.info("archive_month.py completed successfully")
    except Exception as e:
        logger.critical(f"Script failed with exception: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()