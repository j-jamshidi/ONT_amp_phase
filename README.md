# ONT Amplicon Phasing Pipeline

A bioinformatics pipeline for analyzing barcoded amplicon sequences from Oxford Nanopore Technology (ONT) data. The pipeline performs quality control, variant calling, and haplotype phasing to determine whether variants are *cis* (same chromosome) or *trans* (different chromosomes), making it ideal for targeted sequencing studies where specific variants are under investigation.

## Features

- **Quality Control**: Comprehensive read quality assessment including mapping quality (MAPQ≥20), coverage depth statistics, read length distribution (N50), base quality metrics, and alignment identity
- **Variant Calling**: Clair3-based SNP and indel detection with deep neural network accuracy optimized for ONT data
- **Haplotype Phasing**: Dual approach using WhatsHap for robust phasing and HapCUT2 for confirmatory analysis
- **Flexible Input**: Supports single variant QC analysis or dual variant phasing with automatic workflow selection
- **Batch Processing**: Automated analysis of multiple barcoded samples
- **Chimeric Detection**: Identifies and reports reads with discordant phasing patterns indicating potential chimeric events

## Quick Start with Test Data

```bash
# Clone repository
git clone https://github.com/j-jamshidi/ONT_amp_phase.git
cd ONT_amp_phase

# Configure paths in script/run.sh (edit REFERENCE_FASTA, HAPCUT2_PATH)
vim script/run.sh

# Run test dataset
bash script/run.sh test_run

# Check results
ls Runs/test_run/result/barcode*/
```

## Required Folder Structure

The pipeline expects a specific directory organization. Your project should follow this structure:

```
project_root/
├── script/                          # Pipeline scripts and templates
│   ├── run.sh                       # Main pipeline orchestrator
│   ├── basecalling_phasing_amplicon.py  # Two-variant phasing analysis
│   ├── basecalling_QC_amplicon.py   # Single-variant QC analysis
│   └── dummy.vcf                    # VCF template 
└── Runs/                            # Input data directory
    └── <RUN_ID>/                    # Specific run folder (e.g., test_run)
        ├── sample_sheet.csv         # Sample metadata and variant info
        ├── bam_pass/                # ONT basecaller output structure
        │   ├── barcode01/           # Individual barcode directories
        │   │   ├── file1.bam        # BAM files from same barcode
        │   │   └── file2.bam
        │   └── barcode02/
        │       ├── file1.bam
        │       └── file2.bam
        └── result/                  # Output directory (auto-created)
            ├── barcode01/           # Per-sample results
            └── barcode02/
```

## Sample Sheet Format

Create `sample_sheet.csv` in your `Runs/<RUN_ID>/` directory with the following columns:

```csv
Batch,Barcode,Episode,Coordinate,Variant1,Variant2
test,18,NA24385_1,chr4:41637861-41652004,chr4:41638861 C>T,chr4:41651881 G>A
test,19,NA24385_2,chr21:42375008-42389474,chr21:42375588 G>A,chr21:42388818 A>G
test,20,NA24385_3,chr1:100000000-100020000,chr1:100010000 A>G,
```

**Column Descriptions:**
- **Batch**: Run identifier (must match your RUN_ID)
- **Barcode**: Barcode number only (e.g., 18 for barcode18 directory)
- **Episode**: Unique sample identifier (used for output file naming)
- **Coordinate**: Amplicon region in format `chr:start-end` (no commas in coordinates)
- **Variant1**: Primary variant in format `chr:pos REF>ALT` or `chr:pos:REF:ALT` 
- **Variant2**: Secondary variant (leave empty for single-variant QC analysis)

**Analysis Modes:**
- **Two variants**: Performs phasing analysis to determine cis/trans relationship
- **One variant**: Performs QC and variant validation
- **No variants**: Performs amplicon QC and reports all detected variants

## Dependencies

**Required Tools:**
- **samtools** (≥1.10): BAM file manipulation and indexing
- **Clair3** (Docker): Deep learning variant caller for ONT data
- **WhatsHap**: Haplotype phasing using long reads
- **HapCUT2**: Alternative phasing algorithm for confirmation
- **bgzip/tabix**: VCF compression and indexing

**Installation:**
```bash
# Install via conda 
conda install -c bioconda samtools whatshap tabix

# Install HapCUT2 from its GitHub repository
https://github.com/vibansal/HapCUT2

# Install Clair3 Docker
docker pull hkubal/clair3:latest

# Python dependencies
pip install pysam pandas numpy whatshap
```

**System Requirements:**
- Linux/macOS operating system
- Docker (for Clair3)
- Python 3.7+
- Minimum 8GB RAM (16GB recommended for large amplicons)

## Usage

### Initial Setup

1. **Clone repository:**
   ```bash
   git clone https://github.com/j-jamshidi/ONT_amp_phase.git
   cd ONT_amp_phase
   ```

2. **Configure `script/run.sh`:**
   Edit the following variables in `run.sh`:
   ```bash
   REFERENCE_FASTA="/path/to/reference.fa"     # Reference genome (e.g., GRCh38)
   HAPCUT2_PATH="/path/to/hapcut2"             # HapCUT2 installation directory
   CLAIR3_PATH="hkubal/clair3:latest"          # Clair3 Docker image
   ```

3. **Prepare your data:**
   - Organize BAM files in `Runs/<RUN_ID>/bam_pass/barcode##/`
   - Create `sample_sheet.csv` with your sample information

### Running the Pipeline

```bash
# Run with your data
bash script/run.sh <RUN_ID>

# Example with test data
bash script/run.sh test_run

# Monitor progress
tail -f Runs/<RUN_ID>/result/pipeline.log
```

## Output Files

Results are generated in `Runs/<RUN_ID>/result/barcode##/` for each sample:

**Core Output Files:**
- `{Episode}.bam` - Merged amplicon-specific BAM file with index
- `{Episode}.wf_snp.vcf.gz` - Clair3 variant calls (compressed and indexed)
- `{Episode}_report.txt` - Comprehensive analysis report
- `{Episode}_coordinate.bed` - Amplicon region coordinates

**Phasing-Specific Files (two variants only):**
- `{Episode}_phased.bam` - WhatsHap-tagged reads with haplotype information
- `hap2cut_{Episode}` - HapCUT2 phasing results
- `whatshap.log` / `HapCUT2.log` - Tool execution logs

**Report Contents:**

*For QC Analysis (0-1 variants):*
- Amplicon region and length information
- Total reads and QC-passing read counts with pass/fail status
- Depth statistics (median, min, max coverage)
- Read quality metrics (length distribution, N50, MAPQ, base quality)
- Alignment statistics (identity percentage, strand distribution)
- Variant comparison (user-provided vs Clair3-detected variants)

*For Phasing Analysis (2 variants):*
- All QC metrics listed above
- Variant validation results with heterozygosity confirmation
- Spanning read analysis (total, clean, high-quality counts)
- Detailed phasing breakdown (cis vs trans read counts and percentages)
- WhatsHap and HapCUT2 phasing results
- Chimeric read detection and percentage

📖 **[Detailed Report Guide](REPORT_GUIDE.md)** - Complete explanation of report structure and interpretation


## License

MIT License - see LICENSE file for details

## Contact

For questions, issues, or feature requests:
- Email: j.jamshidi@neura.edu.au
- GitHub Issues: https://github.com/j-jamshidi/ONT_amp_phase/issues
