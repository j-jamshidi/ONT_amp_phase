# Amplicon Analysis Pipeline for Oxford Nanopore Sequencing

This repository contains a bioinformatics pipeline designed for analyzing barcoded amplicon sequences generated from Oxford Nanopore Technology (ONT) data. The pipeline is capable of performing quality control, variant calling, and haplotype phasing, specifically tailored for two main purposes:
1. Two-variant phasing to determine whether the variants are on the same chromosome (*cis*) or different chromosomes (*trans*).
2. Single-variant localization and quality control.

## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [How it Works](#how-it-works)
    - [Configuration and Setup](#configuration-and-setup)
    - [Input Data](#input-data)
    - [Pipeline Steps](#pipeline-steps)
- [Output Files](#output-files)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Scripts Overview](#scripts-overview)
- [License](#license)
- [Contact](#contact)

## Introduction

Targeted sequencing of amplicons generated using long-range PCR and sequenced with Oxford Nanopore Technology allows for in-depth analysis of specific genomic regions. This pipeline automates the process of evaluating the quality of these reads, identifying genetic variants within the amplicons (variant calling), and more specifically determining the phase of two pre-determined variants (i.e., whether two variants are on the same chromosome - *cis* - or on different chromosomes - *trans*). The pipeline is flexible, and also supports single-variant localization/QC.

## Features

* **Automated Workflow:** Streamlines the analysis from raw BAMs to a pulished report for phasing or variant localisation. 
* **Quality Control (QC):**
    * Assesses general sequencing quality (read length, mapping quality, base quality, read identity) for each amplicon.
    * Calculates amplicon coverage depth statistics (median, mean, min, max).
    * Filters reads based on user-defined quality thresholds (default QMAP>20) for downstream analysis.
    * Provides a comprehensive QC report.
* **Variant Calling:** Utilizes Clair3 for accurate single nucleotide polymorphism (SNP) and indel calling.
* **Haplotype Phasing:**
    * **WhatsHap:** Employs WhatsHap for robust phasing of heterozygous variants.
    * **HapCUT2:** Incorporates HapCUT2 for an alternative or confirmatory phasing approach.
    * Determines and reports the *cis* or *trans* relationship between two specified variants.
* **Chimeric Read Detection:** Identifies and reports reads that show discordant phasing patterns, indicating potential chimeric events.
* **Flexible Variant Input:** Supports analysis for either one or two specified variants from a sample sheet.
* **Scalability:** Designed to process multiple samples in batches.

## How it Works

The pipeline is orchestrated by a `run.sh` script, which integrates two Python scripts and external bioinformatics tools.

### Configuration and Setup

The `run.sh` script requires initial configuration of base directories, reference genome paths, and tool paths.

* `RUNID`: Identifier for the current sequencing run (the name of the folder with files to analyse / run ID).
* `BASEDIR`: Root directory where raw BAM files and sample sheet are located.
* `WORKDIR`: Directory for storing analysis results and intermediate files.
* `REFERENCE_FASTA`: Path to the reference genome in FASTA format (e.g., GRCh38).
* `CLAIR3_PATH`: Path to the Clair3 variant caller installation. *current version of the script uses clair3 docker. 
* `HAPCUT2_PATH`: Path to the HapCUT2 installation.
* `SCRIPT_PATH`: Path to the `script` directory. The script folder contain the pipeline's scripts and nessesary files. 

### Input Data

The pipeline expects the following input files:

* **Sample Sheet (`sample_sheet.csv`):** A CSV file detailing samples, barcodes, amplicon coordinates, and the variants of interest.
    * `Batch`: Batch number (run ID).
    * `Barcode`: Barcode identifier (e.g., `03`).
    * `Episode`: Unique identifier for the sample (e.g., `NA24385_FY2`).
    * `Coordinate`: Genomic coordinates of the amplicon (e.g., `chr1:236010990-236033398` -must not have comma between numbers).
    * `Variant1`: Details of the first variant (e.g., `chr1:236011853 T>C`or `chr1:236011853:T:C`). If no variant is provided, set to empty.
    * `Variant2`: Details of the second variant (e.g., `chr1:236011853 T>C`or `chr1:236011853:T:C`). If only one variant is provided, set to empty.
* **Raw BAM Files:** Barcoded BAM files organized by barcode within `BASEDIR/bam_pass/`. Each barcode directory should contain BAM files from the same barcode. This is the default structure of the Oxford Nanopore sequening output `BASEDIR/bam_pass/barcode{x}`
* **Dummy VCF Template (`dummy.vcf`):** A template VCF file used to generate sample-specific VCFs based on the `sample_sheet.csv`. This template is located in `SCRIPT_PATH`.

### Pipeline Steps

The `run.sh` script iterates through each sample defined in `sample_sheet.csv` and performs the following operations:

1.  **Prepare Input File (`prepare_input_file`):** Reads the `sample_sheet.csv`, formats barcode numbers (e.g., `barcode01`), converts `Episode` and `EpisodeWES` to uppercase, remove any inconsistencies (e.g extra spaces) and creates a processed `.info` file.
2.  **Create Sample Directory:** A dedicated directory is created for each barcode for results (`WORKDIR/barcodeXX`).
3.  **Prepare VCF File (`prepare_vcf`):**
    * Copies a `dummy.vcf` template to the sample's directory.
    * Populates the template with the `Variant1` and `Variant2` information from the sample sheet.
    * If only one variant is provided, the second variant line is removed from the VCF.
    * The VCF is then sorted.
4.  **Merge BAM Files (`merge_bam_files`):**
    * Merges all BAM files associated with a specific barcode from `BASEDIR/bam_pass/barcodeXX/`.
    * Sorts and indexes the merged BAM.
    * Extracts reads specifically covering the amplicon `Coordinate` into a new BAM file (`${Episode}.bam`).
    * Indexes the amplicon-specific BAM file.
5.  **Create BED File:** Generates a BED file (`${Episode}_coordinate.bed`) containing the amplicon coordinates for use by downstream tools.
6.  **Run Clair3 Variant Calling (`run_clair3`):**
    * Executes Clair3 on the amplicon-specific BAM file and BED file against the reference genome. If needed, Clair3 settings can be adjusted in `run.sh` script.
    * Clair3 performs variant calling and optionally integrates WhatsHap for initial phasing within its output.
    * Intermediate Clair3 directories are removed.
    * The Clair3 output VCF is copied and indexed as `${Episode}.wf_snp.vcf.gz`.
7.  **Run HapCUT2 Phasing (`run_hapcut2`):**
    * **Only if two variants are provided:**
    * Extracts haplotype-informative reads from the BAM file and VCF using `extractHAIRS`.
    * Runs `HAPCUT2` to perform phasing based on the extracted fragments.
8.  **Final Analysis and Cleanup:**
    * **Quality Control (Single Variant / No Variant):** If one or no variants are provided in the sample sheet, the `basecalling_QC_amplicon.py` script is executed. This script performs detailed quality control checks on the amplicon BAM and generates a comprehensive QC report (`${Episode}_report.txt`). It also compares variants found by Clair3 with any user-provided variants.
    * **Phasing Analysis (Two Variants):** If two variants are provided, the `basecalling_phasing_amplicon.py` script is executed. This script performs the following:
        * **Variant Validation:** Checks if the provided variants are adequately covered by reads and exhibit expected heterozygous patterns, flagging insufficient coverage, unexpected alleles, or extreme allele skew. This is to make sure the correct heterozygous variants are provided. 
        * **Quality Control (Phasing Specific):** Filters reads based on mapping quality (MAPQ $\ge 20$) and ensures they span both variant positions (keeps informative reads). It also checks for "clean spanning reads" (reads where both variants can be confidently called as ref or alt).
        * **WhatsHap Phasing:** Runs WhatsHap on the quality-controlled spanning reads to generate a phased VCF and a haplotagged BAM file.
        * **Read-Based Phasing Analysis:** Analyzes the haplotagged BAM file to count reads supporting *cis* (both ref or both alt) and *trans* (one ref, one alt) configurations, providing a percentage breakdown and determining the phase based on read counts.
        * **VCF-Based Phasing Analysis:** Parses the WhatsHap-phased VCF and (if available) HapCUT2-phased VCF to determine the overall phase (Cis/Trans) based on the reported genotypes.
        * Generates a detailed report (`${Episode}_report.txt`) summarizing all QC and phasing results.
    * **Cleanup:** Removes intermediate files and directories.

## Output Files

For each processed sample, the following output files are generated within `WORKDIR/barcodeXX/`:

* **`${Episode}.bam`**: Merged and amplicon-specific BAM file.
* **`${Episode}.bam.bai`**: Index for the amplicon-specific BAM.
* **`${Episode}_coordinate.bed`**: BED file specifying the amplicon region.
* **`${Episode}.wf_snp.vcf.gz`**: GZipped VCF file containing variants called by Clair3.
* **`${Episode}.wf_snp.vcf.gz.tbi`**: Tabix index for the Clair3 VCF.
* **`HapCUT2.log`**: Log file for HapCUT2 execution (if run).
* **`whatshap.log`**: Log file for WhatsHap execution (if run).
* **`${Episode}_report.txt`**: Comprehensive report detailing QC metrics, variant comparison, and phasing results (if applicable).
    * **For QC-only analysis:** Contains amplicon region and length, total/passing QC reads, depth statistics, read length distribution (including N50), mapping quality, base quality, alignment statistics (mean identity), strand distribution, and variant comparison with Clair3 calls.
    * **For Phasing analysis:** Includes variant validation results, detailed QC on spanning reads (total, clean, high-quality spanning reads), WhatsHap and HapCUT2 phasing outcomes, and read-based phasing breakdown (Cis/Trans counts and percentages).
* **`${Episode}.vcf`**: (Temporary, then removed) Dummy VCF used for phasing input.
* **`clean-span-hq.bam`**: (Temporary, then removed) BAM file containing high-quality, clean spanning reads.
* **`${Episode}_phased.bam`**: (Only for two variants) BAM file with reads tagged by WhatsHap for phasing information.
* **`${Episode}_phased.bam.bai`**: Index for the phased BAM.
* **`hap2cut_${Episode}`**: (Only for two variants) HapCUT2 phased VCF (if run).
* **`fragment_${Episode}`**: (Only for two variants) Intermediate file generated by `extractHAIRS`.

## Dependencies

This pipeline relies on several bioinformatics tools and Python libraries.

**External Tools:**

* **`samtools`** (version 1.10 or higher recommended): For BAM file manipulation (merging, sorting, indexing, viewing).
* **`Clair3`**: A deep neural network-based variant caller for ONT data.
* **`WhatsHap`**: A tool for phasing genetic variants using long reads.
* **`HapCUT2`**: A program for constructing haplotypes from sequence data.
* **`bgzip`** and **`tabix`**: For compressing and indexing VCF files.

**Python Libraries:**

* `pysam`: For interacting with BAM and VCF files in Python.
* `pandas`: For data manipulation (used in `basecalling_QC_amplicon.py` for `sample_sheet.csv` processing, though the `run.sh` script currently handles this directly).
* `numpy`: For numerical operations, especially in quality calculations (mean, median, N50).
* `whatshap` (Python library): Underlying library used by WhatsHap tool.

## Usage

1.  **Clone the Repository:**
    ```bash
    git clone https://github.com/j-jamshidi/ONT_amp_phase.git
    cd ONT_amp_phase
    ```
2.  **Configure `run.sh`:** Edit `run.sh` to set `RUNID`, `BASEDIR`, `WORKDIR`, `REFERENCE`, `HAPCUT2_PATH`, and `SCRIPT_PATH` according to your environment.   

3.  **Install Dependencies:** Ensure all external tools (`samtools`, `WhatsHap`, `HapCUT2`, `bgzip`, `tabix`) are installed and accessible. Install Python libraries (`pysam`, `pandas`, `numpy`, `whatshap`) using `pip` or your preferred package manager. Install clair3 docker. 

4.  **Run the Pipeline:**
    ```bash
    bash run.sh <RUN_ID>
    ```
    Replace `<RUN_ID>` with the identifier for your run (e.g., `run.sh my_ont_run`). To run the test set, use `bash run.sh test_run`.


## License

This project is licensed under the MIT License.

## Contact

For questions or issues, please contact j.jamshidi@neura.edu.au.
