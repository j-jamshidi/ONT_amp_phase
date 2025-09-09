# Amplicon Analysis Pipeline for Oxford Nanopore Sequencing

This repository contains a bioinformatics pipeline designed for analyzing barcoded amplicon sequences generated from Oxford Nanopore Technology (ONT) data. The pipeline is capable of performing quality control, variant calling, and haplotype phasing, specifically tailored for two main purposes:
1. Single-variant localization and quality control.
2. Two-variant phasing to determine whether the variants are on the same chromosome (*cis*) or different chromosomes (*trans*).

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

Targeted sequencing of amplicons generated using long-range PCR and sequenced with Oxford Nanopore Technology allows for in-depth analysis of specific genomic regions. This pipeline automates the process of evaluating the quality of these reads, identifying genetic variants within the amplicons, and more specifically determining the phase of two pre-determined variants (i.e., whether two variants are on the same chromosome - *cis* - or on different chromosomes - *trans*). The pipeline is flexible, and also supports single-variant localization/QC.

## Features

* **Automated Workflow:** Streamlines the analysis from raw BAMs to a pulished report for phasing or variant localisation. 
* **Quality Control (QC):**
    * Assesses general sequencing quality (read length, mapping quality, base quality, read identity) for each amplicon.
    * Filters reads based on user-defined quality thresholds (default QMAP>20) for downstream analysis.
    * Provides a comprehensive QC report.
* **Variant Calling:** Utilizes Clair3 for accurate single nucleotide polymorphism (SNP) and indel calling.
* **Haplotype Phasing:**
    * **WhatsHap:** Employs WhatsHap for robust phasing of heterozygous variants.
    * **HapCUT2:** Incorporates HapCUT2 for an alternative or confirmatory phasing approach.
    * Determines and reports the *cis* or *trans* relationship between two specified variants for all reads.
* **Chimeric Read Detection:** Identifies and reports reads that show discordant phasing patterns, indicating potential chimeric events.
* **Flexible Variant Input:** Supports analysis for either one or two specified variants from a sample sheet.

## How it Works

The pipeline is orchestrated by a `run.sh` script, which check/download all necessary tools and integrates all the steps to process amplicon data.   
 

### Input Data
The expected directory structure is as follows:

```
BASEDIR/
├── sample_sheet.csv
└── bam_pass/
   ├── barcode01/
   │   ├── AWC893_pass_1.bam
   │   ├── AWC893_pass_2.bam
   │   └── ...
   ├── barcode02/
   │   ├── AWC894_pass_1.bam
   │   └── ...
   └── ...
```

This structure should be present in your `BASEDIR` before running the pipeline.
* **Sample Sheet (`sample_sheet.csv`):** A CSV file detailing samples, barcodes, amplicon coordinates, and the variants of interest.
    * `Batch`: Batch number (run ID).
    * `Barcode`: Barcode identifier (e.g., `03`).
    * `Episode`: Unique identifier for the sample (e.g., `NA24385_FY2`).
    * `Coordinate`: Genomic coordinates of the amplicon (e.g., `chr1:236010990-236033398` -must not have comma between numbers).
    * `Variant1`: Details of the first variant (e.g., `chr1:236011853 T>C`or `chr1:236011853:T:C`). If no variant is provided, set to empty.
    * `Variant2`: Details of the second variant (e.g., `chr1:236011853 T>C`or `chr1:236011853:T:C`). If only one variant is provided, set to empty.
    * `EpisodeWES`: WES (Whole Exome Sequencing) episode identifier for the sample, set to emty if not available.
* **Raw BAM Files:** Barcoded BAM files organized by barcode within `BASEDIR/bam_pass/`. Each barcode directory should contain BAM files from the same barcode. This is the default structure of the Oxford Nanopore sequening output `BASEDIR/bam_pass/barcode{x}`


### Pipeline Steps

The `run.sh` script iterates through each sample defined in `sample_sheet.csv` and performs the following operations:

1.  **Prepare Input File:** Reads the `sample_sheet.csv` and remove any inconsistencies (e.g extra spaces) and creates a processed `.info` file.
2.  **Create Sample Directory:** A dedicated directory is created for each barcode for results (`WORKDIR/barcodeXX`).
3.  **Prepare VCF File:**
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
    * Executes Clair3 on the amplicon-specific BAM file and BED file against the reference genome.
    * Clair3 performs variant calling and optionally integrates WhatsHap for initial phasing within its output.
    * Intermediate Clair3 directories are removed.
    * The Clair3 output VCF is copied and indexed as `${Episode}.wf_snp.vcf.gz`.    
7.  **Run WhatsHap for Phasing:**
8.  **Run HapCUT2 for Phasing:**
    * **Only if two variants are provided:**
9.  **Final Analysis and Cleanup:**
    * **Quality Control (Single Variant / No Variant):** If one or no variants are provided in the sample sheet, performs detailed quality control checks on the amplicon BAM and generates a comprehensive QC report (`${Episode}_report.txt`). It also compares variants found by Clair3 with any user-provided variants.
    * **Phasing Analysis (Two Variants):** If two variants are provided, performs the following:
        * **Variant Validation:** Checks if the provided variants are adequately covered by reads and exhibit expected heterozygous patterns, flagging insufficient coverage, unexpected alleles, or extreme allele skew.
        * **Quality Control (Phasing Specific):** Filters reads based on mapping quality (MAPQ $\ge 20$) and ensures they span both variant positions. It also checks for "clean spanning reads" (reads where both variants can be confidently called as ref or alt).
        * **Read-Based Phasing Analysis:** Analyzes the haplotagged BAM file to count reads supporting *cis* (both ref or both alt) and *trans* (one ref, one alt) configurations, providing a percentage breakdown and determining the phase based on read counts.
        * **VCF-Based Phasing Analysis:** Parses the WhatsHap-phased VCF and HapCUT2-phased VCF to determine the overall phase (Cis/Trans) based on the reported genotypes.
        * Generates a detailed report (`${Episode}_report.txt`) summarizing all QC and phasing results.
    * **Cleanup:** Removes intermediate files and directories.
    * **Data Upload:** Uploads the processed data for each barcode to an AWS S3 bucket and create an xml file for visualisation. 

## Output Files

For each processed sample, the following output files are generated within `WORKDIR/barcodeXX/`:

* **`${Episode}.bam`**: Merged, coordinate-extracted BAM file for the amplicon.
* **`${Episode}.bam.bai`**: BAM index for `${Episode}.bam`.
* **`${Episode}_coordinate.bed`**: BED file with the amplicon coordinates.
* **`${Episode}.wf_snp.vcf.gz`**: GZipped VCF with Clair3 variant calls.
* **`${Episode}.wf_snp.vcf.gz.tbi`**: Tabix index for the Clair3 VCF.
* **`${Episode}_report.txt`**: Main report with QC metrics, variant comparison, and (if two variants) phasing results.
* **`${Episode}_phased.bam`**: WhatsHap-phased BAM (if two variants).
* **`${Episode}_phased.bam.bai`**: Index for phased BAM.
* **`${Episode}_whathap_phased.vcf.gz`**: WhatsHap-phased VCF (if two variants).
* **`${Episode}_whathap_phased.vcf.gz.tbi`**: Tabix index for WhatsHap VCF.
* **`${Episode}_hapcut2_phased.vcf`**: HapCUT2-phased VCF (if two variants).
* **`hapcut2.log`**: Log file for HapCUT2 (if run).
* **`whatshap.log`**: Log file for WhatsHap (if run).
* **`${Episode}_xml.xml`**: XML file for visualization (always generated).


## Dependencies

This pipeline is **fully containerized** using Docker, requiring minimal local dependencies.

**Docker Images Used:**
* **`hkubal/clair3:latest`**: Official Clair3 container for variant calling.
* **`javadj/ontampip:latest`**: Comprehensive pipeline container containing:


## Usage

1.  **Prerequisites:**
    * Install Docker on your system
    * Ensure Docker daemon is running

2.  **Clone the Repository:**
    * git clone https://github.com/j-jamshidi/ONT_amplicon_phase.git
   cd ONT_amplicon_phase   

3.  **Configure `run.sh`:** 
The `run.sh` script requires initial configuration of base and reference genome paths:

   * `BASEDIR`: Root directory where raw BAM files and sample sheet are located.
   * `WORKDIR`: Directory for storing analysis results and intermediate files.
   * `REFERENCE`: Path to the reference genome in FASTA format (e.g., GRCh38).

3.  **Prepare Input Data:**
    * Place your `sample_sheet.csv` in the `BASEDIR` (e.g.`/EBSDataDrive/ONT/Runs/${RUNID}`).
    * Organize your raw barcoded BAM files in `BASEDIR/bam_pass/barcodeXX/` (The ONT default sequencing output structure).

5.  **Run the Pipeline:**
    ```bash
    bash run.sh <RUN_ID>
    ```
    Replace `<RUN_ID>` with the identifier for your run (e.g., `run.sh my_ont_run`).

**Note:** The pipeline will automatically pull the required Docker images (`hkubal/clair3:latest` and `javadj/ontampip:latest`) on first run.


## License

This project is licensed under the MIT License.

## Container Architecture

The pipeline uses a **hybrid containerization approach** for optimal performance and maintainability:

- **Clair3**: Uses the official `hkubal/clair3:latest` container (maintained by developers)
- **Everything Else**: Uses `javadj/ontampip:latest` containing all other tools and scripts


## Contact

For questions or issues, please open an issue on the GitHub repository.
