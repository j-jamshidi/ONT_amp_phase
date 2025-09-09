# Analysis Report Guide

This guide explains each section of the ONT Amplicon Analysis Pipeline report in detail.

## 📋 Report Header

```
Report for: 25R432204D
------------------------------
Amplicon region: chr2:178637863-178642283
Amplicon length: 4,420 bp
Variant 1: chr2:178640789:C:G
Variant 2: chr2:178639561:A:G
Distance between variants: 1228 bp
```

### Fields Explained
- **Report for**: Sample identifier from the `Episode` column in sample sheet
- **Amplicon region**: Genomic coordinates of the target amplicon (chromosome:start-end)
- **Amplicon length**: Size of the target region in base pairs
- **Variant 1/2**: User-specified variants in format `chromosome:position:reference:alternate`
- **Distance between variants**: Physical distance between the two variants on the genome

---

## ✅ Variant Validation Section

```
=== Variant Validation ===
Variants were validated successfully
Both variants show expected heterozygous patterns in the reads.
```

### What This Means
- **Purpose**: Ensures the specified variants are real and suitable for phasing analysis
- **Validation Criteria**:
  - Sufficient read coverage at variant positions (≥10 reads)
  - Heterozygous pattern (both reference and alternate alleles present)
  - Reasonable allele frequency (typically 20-80% for either allele)
  - No extreme allelic imbalance

### Possible Outcomes
- ✅ **"Variants were validated successfully"**: Both variants pass all checks
- ❌ **"Variant validation failed"**: Issues detected (insufficient coverage, homozygous calls, extreme bias)

---

## 🔍 Quality Control Section

```
=== Quality Control ===
Total reads: 3592
Spanning reads: 2800 (77.95%)
Clean spanning reads: 2511 (89.68% of spanning reads)
High quality spanning reads (MAPQ >= 20): 2509 (99.92% of clean spanning reads)

*QC PASSED (>50 high quality spanning reads)*
```

### Read Categories Explained

#### 1. **Total Reads**
- All reads mapping to the amplicon region
- Includes reads that may not cover both variant positions

#### 2. **Spanning Reads** 
- Reads that cover **both** variant positions
- Essential for phasing analysis
- Percentage shows what fraction of total reads are useful for phasing

#### 3. **Clean Spanning Reads**
- Spanning reads where both variants can be confidently called
- Excludes reads with ambiguous bases (N) or low base quality at variant positions
- High percentage indicates good sequencing quality

#### 4. **High Quality Spanning Reads (MAPQ ≥ 20)**
- Clean spanning reads with high mapping quality
- MAPQ ≥ 20 means <1% chance the read is mismapped
- These reads are used for final phasing analysis

### QC Pass/Fail Criteria
- **PASS**: ≥50 high quality spanning reads
- **FAIL**: <50 high quality spanning reads (insufficient data for reliable phasing)

---

## 🧬 Variant Calling Section

```
===Variant Calling===
Variants called from the amplicon by Clair3:
CHROM	POS	REF	ALT	QUAL	GT
chr2	178638513	G	C	34.3	0/1   
chr2	178638556	CTTTTT	C	7.3	0/1   
chr2	178639561	A	G	53.5	0/1 <- 
chr2	178639839	AACTTAT	A	33.9	0/1   
chr2	178640789	C	G	20.2	0/1 <- 
chr2	178642052	CT	C	6.3	0/1   

Variant Matching Results:
Variant 1: chr2 178639561 A > G - FOUND
Variant 2: chr2 178640789 C > G - FOUND
```

### Clair3 Variant Table
- **CHROM**: Chromosome
- **POS**: Position on chromosome
- **REF**: Reference allele
- **ALT**: Alternate allele
- **QUAL**: Variant quality score (higher = more confident)
- **GT**: Genotype (0/1 = heterozygous, 0/0 = homozygous reference, 1/1 = homozygous alternate)
- **← Arrow**: Indicates user-specified variants found by Clair3

### Variant Matching
- **FOUND**: User variant matches a Clair3-called variant
- **NOT FOUND**: User variant not detected by Clair3 (may indicate low coverage or sequencing issues)

---

## 📊 Phasing Results Section

```
=== Result ===
Total high quality spanning reads: 2509

Detailed categorisation of reads:
Reads with ref allele for both variants (Cis): 894 (35.63%)
Reads with alt allele for both variants (Cis): 1455 (57.99%)
Reads with ref for first, alt for second (Trans): 57 (2.27%)
Reads with alt for first, ref for second (Trans): 103 (4.11%)

Cis reads (both ref or both alt): 2349 (93.62%)
Trans reads (one ref, one alt): 160 (6.38%)

Counting the reads determined the phase as Cis 

Chimeric reads percentage: 6.38% 
```

### Read Categorization

#### **Cis Configuration** (Same Chromosome)
1. **Both Reference (Ref-Ref)**: 
   - Both variants show reference allele on same DNA molecule
2. **Both Alternate (Alt-Alt)**: 
   - Both variants show alternate allele on same DNA molecule

#### **Trans Configuration** (Different Chromosomes)
1. **Ref-Alt**: 
   - First variant = reference, second variant = alternate
2. **Alt-Ref**:
   - First variant = alternate, second variant = reference

### Summary Statistics
- **Total Cis**: variants on same chromosome
- **Total Trans**: variants on different chromosomes
- **Phase Determination**: Majority rule 

### Chimeric Reads
- **Definition**: Reads showing a discordant phase (i.e., not matching the main phase configuration)
- **Interpretation**: 
   - Acceptable percentage (≤40%): Within expected range for technical noise
   - High percentage (>40%): Warning—may indicate technical issues or unreliable phasing results

---

## 🔬 Method Comparison

```
WhatsHap analysis determined phase
HapCUT2 analysis determined phase 
```

### Phasing Methods
- **Read-based counting**: Direct counting of cis vs trans reads (shown above)
- **WhatsHap**: Statistical phasing algorithm using read overlap information
- **HapCUT2**: Alternative phasing algorithm with different statistical approach

- Phasing is considered valid only when confirmed by all three methods


### Quality Flags
- ⚠️ **Low coverage**: <50 spanning reads
- ⚠️ **High noise**: >40% chimeric reads
- ⚠️ **Method disagreement**: Conflicting results between algorithms
- ⚠️ **Variant not found**: User variant not detected by Clair3

