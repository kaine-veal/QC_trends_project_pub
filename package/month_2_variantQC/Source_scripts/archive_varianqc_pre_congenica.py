"""
archive_variantqc_pre_congenica.py

This script processes pre-Congenica Panel sequencing VCF files for a specific month 
and creates a monthly archive CSV.

Pre-Congenica differences:
- VCF pattern: *_filtered_All.vcf (not *haplotypecaller_filtered.vcf.gz)
- Location: May be in different directory structure
- Use this script for months BEFORE Congenica implementation
- Use archive_variant_qc_post_congenica.py for current/post-Congenica data

Monthly archiving workflow:
1. Run this script to archive each pre-Congenica month
2. Run archive_variant_qc_post_congenica.py for post-Congenica months
3. Use plot_variant_trends.py to visualise (works with both)
"""

import gzip
import subprocess
import pathlib
import csv
import re
import argparse
from datetime import datetime

# Set up logging
from utils.logger import create_logger
from utils.config import OUTPUT_ROOT, RESULTS_DIR, apptainer_binds
logger = create_logger('variant_QC.log')

# Archive location for monthly CSVs (same as post-Congenica)
ARCHIVE_BASE = OUTPUT_ROOT / "qc_archive_variantqc"

# Base directory for archived sequencing data
DATA_BASE = RESULTS_DIR


# ============================================================================
# BTCH Date Extraction
# ============================================================================
def extract_btch_date(btch_dir_name):
    """
    Extract run date from BTCH directory name.
    
    Pre-Congenica BTCH directories may have formats like:
    - 250103_BTCH600744_1
    - BTCH600739_1
    
    For dated directories, extract YYMMDD prefix.
    For non-dated, return None (will need to check interop files).
    
    Parameters:
        btch_dir_name (str): BTCH directory name
    
    Returns:
        datetime: Parsed date object, or None if no date in name
    """
    # Try to match YYMMDD prefix (e.g., 250103_BTCH...)
    match = re.match(r'^(\d{6})_BTCH', btch_dir_name)
    if match:
        date_str = match.group(1)
        try:
            year = 2000 + int(date_str[0:2])
            month = int(date_str[2:4])
            day = int(date_str[4:6])
            return datetime(year, month, day)
        except ValueError:
            pass
    
    return None


def extract_date_from_interop(btch_dir):
    """
    Extract run date from interop file in BTCH directory.
    Fallback method when BTCH directory name doesn't include date.
    
    Looks for files like: 250101_MN00000_0001_A000ABCDEF_interop...
    
    Parameters:
        btch_dir (pathlib.Path): Path to BTCH directory
    
    Returns:
        datetime: Parsed date object, or None if not found
    """
    # Look for interop files
    interop_files = list(btch_dir.glob("*_interop*"))
    
    for interop_file in interop_files:
        filename = interop_file.name
        # Extract date from filename (YYMMDD at start)
        date_match = re.match(r'^(\d{6})_', filename)
        if date_match:
            date_str = date_match.group(1)
            try:
                year = 2000 + int(date_str[0:2])
                month = int(date_str[2:4])
                day = int(date_str[4:6])
                return datetime(year, month, day)
            except ValueError:
                continue
    
    return None


# ============================================================================
# bcftools Operations
# ============================================================================

# bcftools is run through Apptainer container
BCFTOOLS_CMD = [
    "apptainer", "exec",
    *apptainer_binds(),
    "docker://quay.io/biocontainers/bcftools:1.21--h3a4d415_1",
    "bcftools"
]


def run_bcftools_command(bcftools_args):
    """
    Execute a bcftools command through Apptainer and return stdout.
    
    Parameters:
        bcftools_args (list): bcftools arguments (e.g., ["stats", "-s", "sample", "file.vcf"])
    
    Returns:
        str: Command stdout
    
    Raises:
        subprocess.CalledProcessError: If command fails
    """
    # Combine base apptainer command with bcftools arguments
    full_command = BCFTOOLS_CMD + bcftools_args
    
    result = subprocess.run(
        full_command,
        check=True,
        text=True,
        capture_output=True
    )
    return result.stdout


def parse_bcftools_stats(stats_output):
    """
    Parse bcftools stats output to extract key metrics.
    
    Extracts:
    - TotalVariants (from SN "number of records")
    - num_snps (for InDel ratio calculation)
    - num_indels (for InDel ratio calculation)
    - TsTv_ratio (from TSTV lines)
    - Het/Hom counts (from PSC lines)
    
    Parameters:
        stats_output (str): Raw bcftools stats output
    
    Returns:
        dict: Parsed metrics
    """
    metrics = {}
    het_count = 0
    hom_alt_count = 0
    tstv_ratio = 0

    for line in stats_output.splitlines():
        # Parse summary statistics (SN lines)
        if line.startswith("SN"):
            fields = line.split("\t")
            if len(fields) < 4:
                continue
                
            label = fields[2].strip().lower()
            value = fields[3].strip()
            
            # Skip missing values
            if value == ".":
                continue
                
            # Convert to numeric
            try:
                numeric_value = float(value) if "." in value else int(value)
            except ValueError:
                continue

            # Extract specific metrics
            if "number of records" in label:
                metrics["TotalVariants"] = numeric_value
            elif "number of snps" in label:
                metrics["num_snps"] = numeric_value
            elif "number of indels" in label:
                metrics["num_indels"] = numeric_value

        # Parse Ts/Tv ratio (TSTV lines)
        elif line.startswith("TSTV"):
            fields = line.split()
            if len(fields) >= 5:
                try:
                    tstv_ratio = float(fields[4])
                except ValueError:
                    pass

        # Parse genotype counts (PSC lines - per-sample counts)
        elif line.startswith("PSC"):
            fields = line.split()
            if len(fields) >= 6:
                try:
                    hom_alt_count = int(fields[4])
                    het_count = int(fields[5])
                except ValueError:
                    pass

    # Calculate Het/Hom ratio
    metrics["HetHom_ratio"] = (
        round(het_count / hom_alt_count, 3) if hom_alt_count > 0 else 0
    )
    
    # Store Ts/Tv ratio
    metrics["TsTv_ratio"] = tstv_ratio
    
    return metrics


def calculate_mean_qual(vcf_path):
    """
    Calculate mean QUAL score using bcftools query.
    
    Parameters:
        vcf_path (pathlib.Path): Path to VCF file
    
    Returns:
        float: Mean QUAL score, or 0 if calculation fails
    """
    try:
        # Query all QUAL values
        qual_output = run_bcftools_command(
            ["query", "-f", "%QUAL\n", str(vcf_path)]
        )
        
        # Parse numeric values (exclude ".")
        qual_values = [
            float(q) for q in qual_output.strip().splitlines()
            if q != "."
        ]
        
        # Calculate mean
        return round(sum(qual_values) / len(qual_values), 2) if qual_values else 0
        
    except Exception:
        return 0


def count_insertions_deletions(vcf_path):
    """
    Count insertions and deletions separately by parsing VCF file.
    
    Parameters:
        vcf_path (pathlib.Path): Path to VCF file
    
    Returns:
        tuple: (num_insertions, num_deletions)
    """
    
    insertions = 0
    deletions = 0
    
    try:
        # Check if file is gzipped
        with open(vcf_path, 'rb') as f:
            magic = f.read(2)
        
        is_gzipped = (magic == b'\x1f\x8b')
        
        # Open file appropriately
        if is_gzipped:
            f = gzip.open(vcf_path, 'rt')
        else:
            f = open(vcf_path, 'r')
        
        for line in f:
            # Skip header lines
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            
            ref = fields[3]
            alt = fields[4]
            
            # Skip multi-allelic sites
            if ',' in alt:
                continue
            
            # Count based on length difference
            if len(alt) > len(ref):
                insertions += 1
            elif len(ref) > len(alt):
                deletions += 1
        
        f.close()
        
    except Exception as e:
        logger.warning(f"Could not count insertions/deletions for {vcf_path.name}: {e}")
        return 0, 0
    
    return insertions, deletions


def get_sample_id(vcf_path):
    """
    Extract sample ID from VCF file.
    
    Parameters:
        vcf_path (pathlib.Path): Path to VCF file
    
    Returns:
        str: Sample ID (W number)
    """
    try:
        sample_output = run_bcftools_command(
            ["query", "-l", str(vcf_path)]
        )
        sample_id = sample_output.strip().splitlines()[0]
        return sample_id
    except subprocess.CalledProcessError:
        # Fallback to filename if sample extraction fails
        return vcf_path.stem


# ============================================================================
# VCF Processing
# ============================================================================
def extract_identifiers_from_path(vcf_path):
    """
    Extract BTCH ID and W folder from VCF file path.
    
    Parameters:
        vcf_path (pathlib.Path): Path to VCF file
    
    Returns:
        tuple: (btch_id, w_folder)
    """
    vcf_string = str(vcf_path)
    
    # Search for BTCH pattern (e.g., BTCH600739_1)
    batch_match = re.search(r"(BTCH\d+_\d+)", vcf_string)
    batch_id = batch_match.group(1) if batch_match else ""
    
    # Search for W folder pattern (e.g., W001)
    w_folder_match = re.search(r"/(W\d+)/", vcf_string)
    w_folder = w_folder_match.group(1) if w_folder_match else ""

    return batch_id, w_folder


def process_vcf(vcf_path):
    """
    Process a single VCF file and extract all QC metrics.
    
    Parameters:
        vcf_path (pathlib.Path): Path to VCF file
    
    Returns:
        dict: Dictionary containing all metrics, or None if processing failed
    """
    # Extract identifiers from path
    batch_id, w_folder = extract_identifiers_from_path(vcf_path)

    # Get sample ID
    sample_id = get_sample_id(vcf_path)

    # Run bcftools stats with sample-specific flag
    try:
        stats_output = run_bcftools_command(
            ["stats", "-s", sample_id, str(vcf_path)]
        )
    except subprocess.CalledProcessError as e:
        logger.error(f"bcftools stats failed for {vcf_path.name}: {e}")
        return None

    # Parse statistics
    metrics = parse_bcftools_stats(stats_output)

    # Calculate mean QUAL
    metrics["MeanQUAL"] = calculate_mean_qual(vcf_path)
    
    # Count insertions and deletions separately for true Ins/Del ratio
    num_insertions, num_deletions = count_insertions_deletions(vcf_path)
    
    # Calculate true Ins/Del ratio (insertions / deletions)
    if num_deletions > 0:
        ins_del_ratio = round(num_insertions / num_deletions, 3)
    else:
        ins_del_ratio = 0

    # Build complete metrics dictionary
    return {
        "BTCH": batch_id,
        "W_folder": w_folder,
        "VCF": vcf_path.name,
        "TotalVariants": metrics.get("TotalVariants", 0),
        "HetHom_ratio": metrics.get("HetHom_ratio", 0),
        "TsTv_ratio": metrics.get("TsTv_ratio", 0),
        "InsDel_ratio": ins_del_ratio,  
        "MeanQUAL": metrics.get("MeanQUAL", 0),
    }


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
    csv_path = month_dir / f"variant_qc_panel_{month}.csv"
    return csv_path.exists()


# ============================================================================
# Main Processing Function
# ============================================================================
def process_month(month: str) -> None:
    """
    Process all pre-Congenica Panel VCF files from a specific month and create archive CSV.
    
    This is the main orchestration function that:
    1. Finds all BTCH directories in the target year
    2. Filters to runs from the specified month (using date from folder name or interop file)
    3. Finds all pre-Congenica Panel VCF files (*_filtered_All.vcf)
    4. Processes each VCF to extract metrics
    5. Writes results to monthly archive CSV
    
    Parameters:
        month (str): Month to process in YYYY-MM format
    """
    logger.info("=" * 60)
    logger.info(f"Processing PRE-CONGENICA Panel variant QC for month: {month}")
    logger.info("=" * 60)

    # Parse target month
    target_year, target_month = map(int, month.split("-"))
    
    # Column names for CSV output (same as post-Congenica for compatibility)
    csv_columns = [
        "BTCH",
        "W_folder",
        "VCF",
        "TotalVariants",
        "HetHom_ratio",
        "TsTv_ratio",
        "InsDel_ratio",
        "MeanQUAL",
    ]
    
    # List to store all results
    all_metrics = []
    
    # Counters for reporting
    total_btch_dirs = 0
    processed_btch = 0
    skipped_btch = 0
    wrong_month_btch = 0
    processed_vcfs = 0
    failed_vcfs = 0
    
    # Find year directory
    year_dir = DATA_BASE / str(target_year)
    if not year_dir.exists():
        logger.error(f"Year directory does not exist: {year_dir}")
        raise FileNotFoundError(f"Year directory not found: {year_dir}")
    
    logger.info(f"Searching for BTCH directories in: {year_dir}")
    
    # Find all BTCH directories in the year
    btch_dirs = sorted([d for d in year_dir.iterdir() if d.is_dir() and "BTCH" in d.name])
    total_btch_dirs = len(btch_dirs)
    logger.info(f"Found {total_btch_dirs} BTCH directories")
    
    # Process each BTCH directory
    for btch_dir in btch_dirs:
        # Try to extract date from directory name first
        run_date = extract_btch_date(btch_dir.name)
        
        # If no date in name, try interop files
        if run_date is None:
            run_date = extract_date_from_interop(btch_dir)
        
        # Skip if we couldn't determine date
        if run_date is None:
            logger.warning(f"Skipping {btch_dir.name} - could not determine run date")
            skipped_btch += 1
            continue
        
        # DEBUG: Log the date we found
        logger.debug(f"{btch_dir.name} -> date found: {run_date.strftime('%Y-%m-%d')}")
        
        # Skip if not target month
        if run_date.year != target_year or run_date.month != target_month:
            logger.info(f"Skipping {btch_dir.name} (date: {run_date.strftime('%Y-%m-%d')}) - not in target month {month}")
            wrong_month_btch += 1
            continue
        
        logger.info(f"Processing BTCH: {btch_dir.name} (run date: {run_date.strftime('%Y-%m-%d')})")
        processed_btch += 1
        
        # PRE-CONGENICA: Find all Panel VCF files with old naming pattern
        # Pattern: W*/*_filtered_All.vcf
        vcf_files = list(btch_dir.glob("W*/*_filtered_All.vcf"))
        
        if not vcf_files:
            logger.warning(f"No pre-Congenica Panel VCF files found in {btch_dir.name}")
            continue
        
        logger.info(f"Found {len(vcf_files)} pre-Congenica Panel VCF files in {btch_dir.name}")
        
        # Process each VCF file
        for vcf_path in vcf_files:
            logger.debug(f"Processing {vcf_path.name}")
            
            try:
                metrics = process_vcf(vcf_path)
                if metrics:
                    all_metrics.append(metrics)
                    processed_vcfs += 1
                else:
                    failed_vcfs += 1
            except Exception as e:
                logger.error(f"Failed to process {vcf_path.name}: {e}")
                failed_vcfs += 1
    
    # Create archive directory if it doesn't exist
    month_dir = ARCHIVE_BASE / month
    month_dir.mkdir(parents=True, exist_ok=True)
    
    # Define output CSV path (same format as post-Congenica for compatibility)
    output_csv = month_dir / f"variant_qc_panel_{month}.csv"
    
    # Write results to CSV
    logger.info(f"Writing {len(all_metrics)} rows to CSV: {output_csv}")
    with open(output_csv, "w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=csv_columns)
        writer.writeheader()
        writer.writerows(all_metrics)
    
    logger.info(f"Metrics successfully saved to {output_csv}")
    
    # Print summary
    logger.info("=" * 60)
    logger.info("Processing Summary:")
    logger.info(f"  - Total BTCH directories found: {total_btch_dirs}")
    logger.info(f"  - Wrong month: {wrong_month_btch}")
    logger.info(f"  - Processed BTCH directories: {processed_btch}")
    logger.info(f"  - Skipped (no date): {skipped_btch}")
    logger.info(f"  - Successfully processed VCFs: {processed_vcfs}")
    logger.info(f"  - Failed VCFs: {failed_vcfs}")
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
            - month: Month to archive (YYYY-MM format)
            - force: Whether to overwrite existing archive
    """
    ap = argparse.ArgumentParser(
        description="Archive PRE-CONGENICA Panel variant QC metrics for a specific month",
    )
    
    ap.add_argument("--month", required=True, help="Month to archive in YYYY-MM format (e.g., 2024-01)")
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
    logger.info("Starting archive_variant_qc_month_pre_congenica.py")
    logger.info("=" * 60)
    
    # Validate month format
    if not validate_month_format(args.month):
        logger.error(f"Invalid month format: {args.month}")
        logger.error("Month must be in YYYY-MM format (e.g., 2024-01)")
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
        csv_path = ARCHIVE_BASE / args.month / f"variant_qc_panel_{args.month}.csv"
        logger.error(f"Existing archive: {csv_path}")
        raise SystemExit("Error: Month already archived (use --force to overwrite)")
    
    # Process the month
    try:
        process_month(args.month)
        logger.info("archive_variant_qc_month_pre_congenica.py completed successfully")
    except Exception as e:
        logger.critical(f"Script failed with exception: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()
