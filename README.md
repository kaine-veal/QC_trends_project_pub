# Program to analyse QC trends for run QC, sample QC and variant QC.
<br>

## Overview

Automated quality control tracking system for Illumina sequencing data with four analysis levels:

1. **Run-Level QC** (Month 0) - Interop metrics from sequencing runs for all runtypes
2. **Sample-Level QC** (Month 1) - FastQC/MultiQC metrics per sample for all runtypes
3. **Sample-Level QC** (Month 1) - Horizontal Coverage Analysis for panels
4. **Variant-Level QC** (Month 2) - VCF metrics for Panels

All systems follow a **two-step workflow**:
1. **Archive** monthly data
2. **Plot** trends from archived data

---

## Prerequisites

**Environment installs**:
```bash
conda env create -f environment.yml
```
```bash
conda activate qc_env
```
```bash
pip install .
```

---

## Configuration

Paths are set through environment variables, so nothing site-specific is hard-coded:

```bash
export QM_RUN_DIR=/path/to/sequencer/runs        # one folder per run
export QM_RESULTS_DIR=/path/to/pipeline/results  # results organised by year
export QM_OUTPUT_ROOT=/path/to/qc_trends         # archives and plots are written here
```

If unset, they default to `data/raw`, `data/results` and `qc_trends` relative to the working directory.

---
<br>
<br>

# Month 0: Run-Level QC (Interop Metrics)
<br>

## Location
Navigate to:
```bash
cd package/month_0_runQC/Source_scripts
```

---

## Archiving

### Usage
```bash
python archive_month.py --rundir $QM_RUN_DIR --month YYYY-MM [--runtype TYPE] [--force]
```

### Arguments
- `--rundir`: Directory containing sequencing runs (usually `$QM_RUN_DIR`)
- `--month`: Month in YYYY-MM format
- `--runtype`: *Optional* - Filter by `RACP`, `Panel`, or `Genotyping` (faster processing)
- `--force`: Overwrite existing archive

### Examples
```bash
# Archive all runtypes for January 2025
python archive_month.py --rundir $QM_RUN_DIR --month 2025-01

# Archive only Panel runs (faster)
python archive_month.py --rundir $QM_RUN_DIR --month 2025-01 --runtype Panel

# Overwrite existing archive
python archive_month.py --rundir $QM_RUN_DIR --month 2025-01 --force
```

### Output
Archives saved to:
```bash
$QM_OUTPUT_ROOT/qc_archive_runqc/YYYY-MM/
```

---

## Plotting

### Usage
```bash
python plot_trends.py --months YYYY-MM[,YYYY-MM...] --metrics METRIC [--runtype TYPE] --out FILENAME
```

### Arguments
- `--months`: Comma-separated months (e.g., `2025-01,2025-02`)
- `--metrics`: Metrics to plot (see table below)
- `--runtype`: *Optional for runtime/yield* - `RACP`, `Panel`, or `Genotyping`
- `--out`: Output filename
- `--list-months`: List available archived months

### Available Metrics

| Metric | Aliases | Runtype Required | Description |
|--------|---------|------------------|-------------|
| `Q30_%` | `q30` | Yes (single only) | % bases ≥Q30 |
| `PhiX_Error_Rate_%` | `error_rate`, `phix_error` | Yes (single only) | PhiX alignment error rate |
| `Total_Runtime_hours` | `runtime`, `runtime_hours` | No | Run duration |
| `Total_Yield_Gb` | `yield` | No | Total sequencing yield |

### Examples
```bash
# List available months
python plot_trends.py --list-months

# Plot Q30 for RACP (single month)
python plot_trends.py --months 2025-01 --metrics Q30 --runtype RACP --out racp_jan.png

# Plot multiple metrics for Panel (Q1 2025)
python plot_trends.py --months 2025-01,2025-02,2025-03 --metrics Q30 error_rate --runtype Panel --out panel_q1.png

# Plot universal metrics for all runtypes
python plot_trends.py --months 2025-01,2025-02 --metrics runtime yield --out all_jan_feb.png
```

### Output
Plots saved to:
```bash
$QM_OUTPUT_ROOT/month_0_output/plot_output/
```

---
<br>
<br>

# Month 1: Sample-Level QC (FastQC/MultiQC)
<br>

## Location
```bash
cd package/month_1_sampleQC/Source_scripts
```

---

## Step 1: Archiving: Generate FastQC/MultiQC Data (Monthly Requirement)

### 1.1 Create Monthly Directory

**For Genotyping**:
```bash
cd $QM_OUTPUT_ROOT/qc_archive/genotyping
```
```
mkdir 2025-12
```
```
cd 2025-12
```
```
mkdir fastqc_output
```
```
mkdir multiqc_output
```

**For RACP/Panels**: Follow same pattern in respective directories

### 1.2 Run FastQC Script

Navigate to scripts:
```bash
cd <fastqc_scripts_dir>
```

Edit script (example for Genotyping):
```bash
nano genotyping_fastqc_by_run.sh
```

**Required edits**:
- **Line 8**: Update output path
  ```bash
  OUTPUT_BASE="$QM_OUTPUT_ROOT/qc_archive/genotyping/2025-12/fastqc_output"
  ```
- **Lines 11-12**: Set date range for the month
- **Line 135**: Update MultiQC output
  ```bash
  MULTIQC_OUT="$QM_OUTPUT_ROOT/qc_archive/genotyping/2025-12/multiqc_output"
  ```

Run:
```bash
./genotyping_fastqc_by_run.sh
```

---

## Step 2: Plotting: Generate Reports

### Usage
```bash
python extract_and_plot.py --run-type TYPE --months YYYY-MM[,YYYY-MM...] [--metrics METRIC1 METRIC2...]
```

### Arguments
- `--run-type`: `racp`, `genotyping`, or `panels`
- `--months`: Comma-separated months
- `--metrics`: *Optional* - Filter specific metrics (see table below)
- `--list-months`: List available archived months

### Available Metrics

Each metric has `_fail_proportion` and `_warn_proportion` variants:

| Metric Base Name | Description |
|------------------|-------------|
| `per_base_sequence_quality` | Base quality across read positions |
| `per_tile_sequence_quality` | Quality variation between flow cell tiles |
| `per_sequence_quality_scores` | Overall read quality distribution |
| `per_base_n_content` | Ambiguous base calls (N) |
| `adapter_content` | Adapter contamination |
| `overrepresented_sequences` | Duplicate/overrepresented sequences |

**⚠️ Note**: Must specify full name with `_fail_proportion` or `_warn_proportion` suffix

### Examples
```bash
# List available months
python extract_and_plot.py --run-type racp --list-months

# Single month (all metrics)
python extract_and_plot.py --run-type racp --months 2025-01

# Multiple months (trend analysis)
python extract_and_plot.py --run-type racp --months 2025-01,2025-02,2025-03

# Specific metrics
python extract_and_plot.py --run-type genotyping --months 2025-01 --metrics per_base_sequence_quality_fail_proportion adapter_content_warn_proportion
```

### Output
Files saved to:
```
$QM_OUTPUT_ROOT/month_1_output/{runtype}_output/
```

**Files generated**:
- `{runtype}_qc_trends_{month}.csv` - All metrics data
- `{runtype}_qc_trends_{month}.png` - Visualization
- `{runtype}_qc_trends_{month}_table.xlsx` - Excel summary

---
<br>
<br>

# Month 1: Horizontal Coverage Analysis
<br>

### Location
Navigate to:
```bash
cd package/month_1_sampleQC/Source_scripts
```

### Usage
```bash
python horcov_analysis.py [--months YYYY-MM,YYYY-MM...]
```

### Arguments
- `--months`: *Optional* - Specific months to analyse (e.g., `2025-01,2025-05`)
- *No arguments*: Analyses all data from current year

### Examples
```bash
# Analyse current year
python horcov_extraction_and_plot.py

# Specific months
python horcov_extraction_and_plot.py --months "2025-01,2025-05,2025-12"
```

### What It Does
- Analyses exon coverage at ≥20× depth
- Identifies samples below 90% coverage threshold
- Generates per-run boxplots with sample overlays
- Shows median and mean coverage per run

### Output
Files saved to:
```
$QM_OUTPUT_ROOT/month_1_output_HorCov/
```

**Files generated**:
- `coverage_all_samples_[months]_TIMESTAMP.csv` - All coverage data
- `outliers_below_90pct_[months]_TIMESTAMP.csv` - Samples <90% coverage
- `coverage_20x_[months]_TIMESTAMP.png` - Boxplot visualization

---
<br>
<br>

# Month 2: Variant-Level QC (Panel Sequencing)
<br>

## Location
```bash
cd package/month_2_variantQC/Source_scripts
```

---

## Quick Start

**Archive pre-Congenica** (April-September 2025):
```bash
python archive_variantqc_pre_congenica.py --month 2025-04
```

**Archive post-Congenica** (October 2025 onwards):
```bash
python archive_variantqc_post_congenica.py --month 2025-10
```

**Plot trends**:
```bash
python plot_variantqc.py --months 2025-04,2025-10 --metric InsDel_ratio --out insdel.png
```

---

## Timeline

| Period | VCF Format | Script to Use |
|--------|-----------|---------------|
| **2025-01 to 2025-09** | `*_filtered_All.vcf` | Pre-Congenica |
| **2025-10 onwards** | `*haplotypecaller_filtered.vcf.gz` | Post-Congenica |

---

## Archiving

### Pre-Congenica Script

**Usage**:
```bash
python archive_variantqc_pre_congenica.py --month YYYY-MM [--force]
```

**Example**:
```bash
# Single month
python archive_variantqc_pre_congenica.py --month 2025-04
```

### Post-Congenica Script

**Usage**:
```bash
python archive_variantqc_post_congenica.py --month YYYY-MM [--force]
```

**Example**:
```bash
# Single month
python archive_variantqc_post_congenica.py --month 2025-10
```

### Output
Archives saved to:
```
$QM_OUTPUT_ROOT/variant_qc_archive/YYYY-MM/
```

---

## Plotting

### Usage
```bash
python plot_variantqc.py --months YYYY-MM[,YYYY-MM...] --metric METRIC --out FILENAME
```

### Arguments
- `--months`: Comma-separated months
- `--metric`: Metric to plot (see table below)
- `--out`: Output filename
- `--list-months`: List available archived months

### Available Metrics

| Metric | Typical Range | Description |
|--------|---------------|-------------|
| `TotalVariants` | 200-500 | Total variants per sample |
| `HetHom_ratio` | 1.5-2.5 | Heterozygous/Homozygous ratio |
| `TsTv_ratio` | 2.0-2.5 | Transition/Transversion ratio |
| `InsDel_ratio` | 0.5-2.0 | Insertions/Deletions ratio  |
| `MeanQUAL` | 6k-14k | Average Phred quality score |


### Examples
```bash
# List available months
python plot_variantqc.py --list-months

# Single metric, multiple months
python plot_variantqc.py --months 2025-04,2025-05,2025-10 --metric InsDel_ratio --out insdel_trend.png
```

### Output
Plots saved to:
```
$QM_OUTPUT_ROOT/variant_qc_output/plots/
```

---



### Why MeanQUAL is So High
HaplotypeCaller uses pair-HMM model generating strong confidence scores. Values of 6,000-14,000 are **normal and expected**.

### InsDel Ratio Calculation
Parses VCF directly to count insertions vs deletions:
```python
if len(alt) > len(ref): insertions += 1
elif len(ref) > len(alt): deletions += 1
InsDel_ratio = insertions / deletions
```

---
<br>
<br>

# Known Issues & To-Do
<br>

## Month 0 (Run QC)
- Unknown runtypes (no SampleSheet) are logged but still plotted - should be excluded

## Month 1 (Sample QC)
- Gaps in FASTQ retention leave gaps in the sample QC trends

## Month 2 (Variant QC)
- Add dbSNP annotation feature
- Genotyping variant QC not implemented (amplicon-specific grouping needed)
- Consider log₁₀ scale for MeanQUAL plots

## All Systems
- Make interactive plots (e.g., using ChronQC)
- Auto-size graphs to number of runs being visualized
- Consolidate monthly maintenance into single command across all systems

---
