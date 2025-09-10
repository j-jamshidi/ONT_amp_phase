# ONT Amplicon Phasing Pipeline

This repository contains a bioinformatics pipeline designed for analyzing barcoded amplicon sequences generated from Oxford Nanopore Technology (ONT) data. The pipeline is capable of performing quality control, variant calling, and haplotype phasing, specifically tailored for two main purposes:
1. Single-variant localization and quality control.
2. Two-variant phasing to determine whether the variants are on the same chromosome (*cis*) or different chromosomes (*trans*).

## Table of Contents
- [Introduction](#introduction)
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

Targeted sequencing of amplicons generated using long-range PCR and sequenced with Oxford Nanopore Technology allows for in-depth analysis of specific genomic regions. This pipeline automates the process of evaluating the quality of these reads, identifying genetic variants within the amplicons, and more specifically determining the phase of two pre-determined variants (i.e., whether two variants are on the same chromosome - *cis* - or on different chromosomes - *trans*). The pipeline is flexible, and also supports single-variant localization/QC.


## How it Works

The pipeline is orchestrated by the `ontphase` script, which checks/downloads all necessary tools and integrates all the steps to process amplicon data.

### Input Data
The expected directory structure is as follows:

```
INPUT_DIR/
├── sample_sheet.csv
└── bam_pass/
   ├── barcode01/
   │   ├── sample_pass_1.bam
   │   ├── sample_pass_2.bam
   │   └── ...
   ├── barcode02/
   │   ├── sample_pass_1.bam
   │   └── ...
   └── ...
```

* **Sample Sheet (`sample_sheet.csv`):** A CSV file detailing samples, barcodes, amplicon coordinates, and the variants of interest.
    * `Batch`: Batch number (run ID).
    * `Barcode`: Barcode identifier (e.g., `03`).
    * `Episode`: Unique identifier for the sample (e.g., `NA24385_FY2`).
    * `Coordinate`: Genomic coordinates of the amplicon (e.g., `chr1:236010990-236033398` - must not have comma between numbers).
    * `Variant1`: Details of the first variant (e.g., `chr1:236011853 T>C` or `chr1:236011853:T:C` - must not have comma between numbers). If no variant is provided, set to empty.
    * `Variant2`: Details of the second variant (e.g., `chr1:236011853 T>C` or `chr1:236011853:T:C` - must not have comma between numbers). If only one variant is provided, set to empty.
* **Raw BAM Files:** Barcoded BAM files organized by barcode within `INPUT_DIR/bam_pass/`. Each barcode directory should contain BAM files from the same barcode. This is the default structure of Oxford Nanopore sequencing output.


### Pipeline Steps

The `ontphase` script iterates through each sample defined in `sample_sheet.csv` and performs the following operations:

1.  **Prepare Input File:** Reads the `sample_sheet.csv`, removes any inconsistencies (e.g., extra spaces), and creates a processed `.info` file.
2.  **Create Sample Directory:** A dedicated directory is created for each barcode for results (`OUTPUT_DIR/barcodeXX`).
3.  **Prepare VCF File:**
    * Copies a `dummy.vcf` template to the sample's directory.
    * Populates the template with the `Variant1` and `Variant2` information from the sample sheet.
    * If only one variant is provided, the second variant line is removed from the VCF.
    * The VCF is then sorted.
4.  **Merge BAM Files:**
    * Merges all BAM files associated with a specific barcode from `INPUT_DIR/bam_pass/barcodeXX/`.
    * Sorts and indexes the merged BAM.
    * Extracts reads specifically covering the amplicon `Coordinate` into a new BAM file (`${Episode}.bam`).
    * Indexes the amplicon-specific BAM file.
5.  **Create BED File:** Generates a BED file (`${Episode}_coordinate.bed`) containing the amplicon coordinates for use by downstream tools.
6.  **Run Clair3 Variant Calling:**
    * Executes Clair3 on the amplicon-specific BAM file and BED file against the reference genome.
    * Clair3 performs variant calling.
    * Intermediate Clair3 directories are removed.
    * The Clair3 output VCF is copied and indexed as `${Episode}.wf_snp.vcf.gz`.
7.  **Run WhatsHap for Phasing:** (Only if two variants are provided)
8.  **Run HapCUT2 for Phasing:** (Only if two variants are provided)
9.  **Final Analysis and Cleanup:**
    * **Quality Control (Single Variant / No Variant):** If one or no variants are provided in the sample sheet, performs detailed quality control checks on the amplicon BAM and generates a comprehensive QC report (`${Episode}_report.txt`). It also compares variants found by Clair3 with any user-provided variants.
    * **Phasing Analysis (Two Variants):** If two variants are provided, performs the following:
        * **Variant Validation:** Checks if the provided variants are adequately covered by reads and exhibit expected heterozygous patterns, flagging insufficient coverage, unexpected alleles, or extreme allele skew.
        * **Quality Control (Phasing Specific):** Filters reads based on mapping quality (MAPQ ≥ 20) and ensures they span both variant positions. It also checks for "clean spanning reads" (reads where both variants can be confidently called as ref or alt).
        * **Read-Based Phasing Analysis:** Analyzes the haplotagged BAM file to count reads supporting *cis* (both ref or both alt) and *trans* (one ref, one alt) configurations, providing a percentage breakdown and determining the phase based on read counts.
        * **VCF-Based Phasing Analysis:** Parses the WhatsHap-phased VCF and HapCUT2-phased VCF to determine the overall phase (Cis/Trans) based on the reported genotypes.
        * Generates a detailed report (`${Episode}_report.txt`) summarizing all QC and phasing results.
    * **Cleanup:** Removes intermediate files and directories. 

## Output Files

For each processed sample, the following output files are generated within `OUTPUT_DIR/barcodeXX/`:

* **`${Episode}.bam`**: Merged, coordinate-extracted BAM file for the amplicon.
* **`${Episode}.bam.bai`**: BAM index for `${Episode}.bam`.
* **`${Episode}_coordinate.bed`**: BED file with the amplicon coordinates.
* **`${Episode}.wf_snp.vcf.gz`**: GZipped VCF with Clair3 variant calls.
* **`${Episode}.wf_snp.vcf.gz.tbi`**: Tabix index for the Clair3 VCF.
* **`${Episode}_report.txt`**: Main report with QC metrics, variant comparison, and (if two variants) phasing results. See [Report Guide](REPORT_GUIDE.md) for detailed explanation of report contents.
* **`${Episode}_phased.bam`**: WhatsHap-phased BAM (if two variants).
* **`${Episode}_phased.bam.bai`**: Index for phased BAM.
* **`${Episode}_Phased.vcf.gz`**: WhatsHap-phased VCF (if two variants).
* **`${Episode}_Phased.vcf.gz.tbi`**: Tabix index for WhatsHap VCF.
* **`hap2cut_${Episode}`**: HapCUT2-phased output (if two variants).
* **`HapCUT2.log`**: Log file for HapCUT2 (if run).
* **`whatshap.log`**: Log file for WhatsHap (if run).
* **`pipeline.log`**: Detailed pipeline execution log for each sample.


## Dependencies

This pipeline is **fully containerized** using Docker, requiring minimal local dependencies.

**Platform Compatibility:**
* Tested on Linux and macOS (including Apple Silicon M-series chips)
* Requires Docker Desktop or Docker Engine

**Docker Images Used:**
* **`hkubal/clair3:latest`**: Official Clair3 container for variant calling.
* **`javadj/ontampip:latest`**: Comprehensive pipeline container containing:
  * WhatsHap for haplotype phasing
  * HapCUT2 for alternative phasing
  * Python scripts for quality control and analysis
  * All necessary bioinformatics tools


## Usage

1.  **Prerequisites:**
    * Install Docker on your system
    * Ensure Docker daemon is running

2.  **Clone the Repository:**
    ```bash
    git clone https://github.com/j-jamshidi/ONT_amplicon_phase.git
    cd ONT_amplicon_phase
    ```

3.  **Make the script executable:**
    ```bash
    chmod +x ontphase
    ```

4.  **Prepare Input Data:**
    * Place your `sample_sheet.csv` in your input directory.
    * Organize your raw barcoded BAM files in `INPUT_DIR/bam_pass/barcodeXX/` (The ONT default sequencing output structure).

5.  **Run the Pipeline:**
    ```bash
    ./ontphase -i <INPUT_DIR> -o <OUTPUT_DIR> -r <REFERENCE_GENOME>
    ```
    
    **Parameters:**
    * `-i, --input`: Input directory containing sample_sheet.csv and bam_pass/ subdirectory
    * `-o, --output`: Output directory for results
    * `-r, --ref`: Path to reference genome FASTA file (e.g., GRCh38)
    * `-h, --help`: Show help message
    
    **Example:**
    ```bash
    ./ontphase -i /path/to/run_data -o /path/to/results -r /path/to/GRCh38.fna
    ```

**Note:** The pipeline will automatically pull the required Docker images (`hkubal/clair3:latest` and `javadj/ontampip:latest`) on first run.

## Test Run

A test dataset is included in the repository to verify the pipeline installation and functionality:

1. **After cloning the repository and making the script executable:**
   ```bash
   ./ontphase -i test_run -o test_results -r /path/to/your/reference.fna
   ```

2. **The test dataset includes:**
   - Sample sheet with two test samples (barcode 18 with one variant for localisation and barcode19 with two varaints for phasing analysis)
   - Barcoded BAM files in the expected directory structure
  
3. **Expected output:**
   - Results will be generated in `test_results/` directory
   - Each barcode will have its own subdirectory with analysis results
   - Check the `*_report.txt` files for QC metrics and/or phasing results

**Note:** You will need to provide your own reference genome file (e.g., GRCh38.fna) as this is not included in the repository due to file size constraints.

## License

This project is licensed under the MIT License.

## Contact

For questions or issues, please contact j.jamshidi@neura.edu.au.
