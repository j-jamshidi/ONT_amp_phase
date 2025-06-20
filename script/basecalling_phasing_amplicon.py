import pysam
import argparse
import subprocess
from collections import Counter
import whatshap
from whatshap.vcf import VcfReader, VariantTable
from whatshap.align import edit_distance
from whatshap.variants import ReadSetReader
import os
import sys
from typing import Optional, Tuple, List
from pathlib import Path

#Change the path to your own reference file
REFERENCE_FASTA = "/Users/javadjamshidi/Desktop/Refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

def compare_variants(dummy_vcf: str, variant_calling_vcf: str) -> List[Tuple[str, int, str, str]]:

    # Read dummy VCF variants
    dummy_variants = []
    with open(dummy_vcf, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    chrom, pos, _, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
                    dummy_variants.append((chrom, pos, ref, alt))
                else:
                    print(f"Skipping malformed line in dummy VCF: {line.strip()}")

    # Read variant calling VCF variants
    variant_calling_variants = []
    try:
        vcf = pysam.VariantFile(variant_calling_vcf)
        for record in vcf.fetch():
            # Only add variants with valid alternate alleles (not "." or None)
            if record.alts and record.alts[0] and record.alts[0] != ".":
                variant_calling_variants.append((record.chrom, record.pos, record.ref, record.alts[0]))
        
        print(f"Number of variants called from the amplicon: {len(variant_calling_variants)}")
    except Exception as e:
        print(f"Error reading variant calling VCF: {e}")
        return []

    # Compare variants
    matched_variants = []
    for dummy_var in dummy_variants:
        found = False
        for vc_var in variant_calling_variants:
            if (dummy_var[0] == vc_var[0] and  # chromosome
                dummy_var[1] == vc_var[1] and   # position
                dummy_var[2] == vc_var[2] and   # reference allele
                dummy_var[3] == vc_var[3]):     # alternate allele
                found = True
                matched_variants.append(dummy_var)
                break
        
        if not found:
            matched_variants.append(None)
    
    return matched_variants

def write_variant_comparison_results(dummy_vcf: str, variant_calling_vcf: str, output_file: str):
    # Ensure dummy VCF exists
    if not os.path.exists(dummy_vcf):
        print(f"Error: Dummy VCF file does not exist at {dummy_vcf}")
        return
    
    # Ensure variant calling VCF exists
    if not os.path.exists(variant_calling_vcf):
        print(f"Error: Variant calling VCF file does not exist at {variant_calling_vcf}")
        return

    try:
        matched_variants = compare_variants(dummy_vcf, variant_calling_vcf)
        
        with open(output_file, 'a') as f:  # Changed to append mode

            f.write("Variants to be phased:\n")
            # Read dummy VCF variants and store them
            variants = []
            dummy_vcf_content = []
            with open(dummy_vcf, 'r') as dummy_f:
                for line in dummy_f:
                    if not line.startswith('#'):
                        parts = line.strip().split('\t')
                        if len(parts) >= 5:
                            chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
                            variants.append((chrom, int(pos), ref, alt))
                            dummy_vcf_content.append(f"{chrom}\t{pos}\t{ref}\t{alt}\n")

            # Create a set of them for matching
            dummy_variants_set = set()
            with open(dummy_vcf, 'r') as dummy_f:
                for line in dummy_f:
                    if not line.startswith('#'):
                        parts = line.strip().split('\t')
                        if len(parts) >= 10:
                            chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
                            dummy_variants_set.add((chrom, int(pos), ref, alt))
            
            # Write the variants and distance information
            for line in dummy_vcf_content:
                f.write(line)
            
            if len(variants) >= 2:
                f.write(f"Distance between the variants is {variants[1][1] - variants[0][1]:,} bp\n")

            f.write("\n===Variant Calling===\n")           
            f.write("Variants called from the amplicon by Clair3:\n")
            f.write("CHROM\tPOS\tREF\tALT\tQUAL\tGT\n")
            # Read variant calling VCF variants
            vcf = pysam.VariantFile(variant_calling_vcf)
            for record in vcf.fetch():
                if record.alts and record.alts[0] and record.alts[0] != ".":
                    qual = f"{record.qual:.1f}" if record.qual is not None else "."
                    gt = record.samples[0]['GT']
                    gt_str = "/".join(str(x) for x in gt) if gt is not None else "."
                    
                    # Check if this variant matches any in the dummy VCF
                    is_matched = (record.chrom, record.pos, record.ref, record.alts[0]) in dummy_variants_set
                    
                    # Add ">" for matched variants
                    indicator = "<-" if is_matched else " "
                    f.write(f"{record.chrom}\t{record.pos}\t{record.ref}\t{record.alts[0]}\t{qual}\t{gt_str} {indicator} \n")
            
            f.write("\nVariant Matching Results:\n")
            for i, variant in enumerate(matched_variants, 1):
                if variant:
                    f.write(f"Variant {i}: {variant[0]} {variant[1]} {variant[2]} > {variant[3]} - FOUND\n")
                else:
                    f.write(f"Variant {i}: NOT FOUND in variant calling VCF\n")
            
            # Summary statistics
            found_count = sum(1 for v in matched_variants if v)
            not_found_count = sum(1 for v in matched_variants if v is None)
    
    except Exception as e:
        print(f"Error during variant comparison: {e}")
        import traceback
        traceback.print_exc()

def normalize_allele(allele: str) -> str:
    """Normalize allele to uppercase and replace empty alleles with <DEL>."""
    if allele == "":
        return "<DEL>"
    return allele.upper()

def get_allele(read, variant) -> Optional[str]:
    """
    Determine which allele is present in a read at a given variant position.
    
    Args:
        read: pysam.AlignedSegment object
        variant: pysam.VariantRecord object
    
    Returns:
        'ref', 'alt', or None if the allele cannot be determined
    """
    pos = variant.pos - 1  # Convert to 0-based position
    
    # Check if the variant position is covered by the read
    if pos < read.reference_start or pos >= read.reference_end:
        return None
    
    # Get the reference and alternate alleles
    ref_allele = variant.ref.upper()
    alt_allele = variant.alts[0].upper() if variant.alts else ""
    
    # Get the aligned pairs with sequence
    aligned_pairs = read.get_aligned_pairs(with_seq=True)
    
    # Extract the variant region with context
    context_size = max(len(ref_allele), len(alt_allele)) + 10  # Increased context size
    variant_region = []
    ref_pos_idx = None
    
    # Find the variant position and collect context
    for idx, (qpos, rpos, base) in enumerate(aligned_pairs):
        if rpos == pos:
            ref_pos_idx = idx
            # Collect preceding context
            start_idx = max(0, idx - context_size//2)
            variant_region.extend(aligned_pairs[start_idx:idx])
            break
    
    if ref_pos_idx is None:
        return None
    
    # Collect following context
    end_idx = min(len(aligned_pairs), ref_pos_idx + context_size//2)
    variant_region.extend(aligned_pairs[ref_pos_idx:end_idx])
    
    # Extract sequences for comparison
    read_seq = ""
    ref_seq = ""
    read_positions = []
    ref_positions = []
    
    for qpos, rpos, base in variant_region:
        if qpos is not None and rpos is not None:
            read_seq += read.query_sequence[qpos]
            ref_seq += base.upper() if base else 'N'
            read_positions.append(qpos)
            ref_positions.append(rpos)
        elif qpos is None and rpos is not None:
            # Deletion in read
            read_seq += '-'
            ref_seq += base.upper() if base else 'N'
            read_positions.append(None)
            ref_positions.append(rpos)
        elif qpos is not None and rpos is None:
            # Insertion in read
            read_seq += read.query_sequence[qpos]
            ref_seq += '-'
            read_positions.append(qpos)
            ref_positions.append(None)
    
    # Handle different variant types
    if len(ref_allele) != len(alt_allele):
        # Indel case
        variant_start_idx = ref_positions.index(pos) if pos in ref_positions else -1
        if variant_start_idx == -1:
            return None
            
        # For insertions 
        if len(alt_allele) > len(ref_allele):
            # Extract the sequence around the insertion point
            insertion_length = len(alt_allele) - len(ref_allele)
            sequence_window = 5  # Consider surrounding context
            
            # Get the read sequence at and after the variant position
            read_variant_seq = ""
            current_idx = variant_start_idx
            insertion_bases = 0
            context_bases = 0
            total_bases_needed = insertion_length + sequence_window
            
            while current_idx < len(read_positions) and len(read_variant_seq) < total_bases_needed:
                if read_positions[current_idx] is not None:
                    read_variant_seq += read_seq[current_idx]
                    
                    # Count insertion bases (those without corresponding reference position)
                    if ref_positions[current_idx] is None:
                        insertion_bases += 1
                    else:
                        context_bases += 1
                current_idx += 1
            
            # Check if we have enough bases for comparison
            if len(read_variant_seq) >= insertion_length:
                # Calculate alignment scores against both alleles
                alt_score = 0
                ref_score = 0
                
                # Compare with alt allele
                for i in range(min(len(read_variant_seq), len(alt_allele))):
                    if read_variant_seq[i] == alt_allele[i]:
                        alt_score += 1
                
                # Compare with ref allele
                for i in range(min(len(read_variant_seq), len(ref_allele))):
                    if read_variant_seq[i] == ref_allele[i]:
                        ref_score += 1
                
                # Add weight to the insertion detection
                if insertion_bases >= insertion_length * 0.8:  # Allow for some sequencing errors
                    alt_score += 2  # Bonus for having the right number of inserted bases
                
                # Calculate match percentages
                alt_match_percent = alt_score / len(alt_allele)
                ref_match_percent = ref_score / len(ref_allele)
                
                # Make the final call
                if alt_match_percent >= 0.8 and insertion_bases >= insertion_length * 0.8:
                    return 'alt'
                elif ref_match_percent >= 0.8 and insertion_bases == 0:
                    return 'ref'
                
        # For deletions 
        else:
            deletion_length = len(ref_allele) - len(alt_allele)
            expected_ref_positions = list(range(pos, pos + len(ref_allele)))
            
            # Count missing bases in the read
            missing_bases = 0
            for exp_pos in expected_ref_positions:
                if exp_pos not in ref_positions or read_positions[ref_positions.index(exp_pos)] is None:
                    missing_bases += 1
            
            # Calculate match scores for both alleles
            ref_match_score = len(ref_allele) - missing_bases
            alt_match_score = missing_bases
            
            # Add additional check for partial deletions
            if missing_bases >= deletion_length * 0.8:  # Allow some mismatches
                return 'alt'
            elif ref_match_score >= len(ref_allele) * 0.8:  # Allow some mismatches
                return 'ref'
    
    # Handle SNVs 
    else:
        query_pos = aligned_pairs[ref_pos_idx][0]
        if query_pos is None:
            return None
        
        read_base = read.query_sequence[query_pos]
        if read_base == ref_allele:
            return 'ref'
        elif read_base == alt_allele:
            return 'alt'
    
    return None
    

def run_whatshap(bam_file, vcf_file, output_bam, reference_fasta):
    """
    Run WhatsHap for phasing variants and tagging reads.
    Ensures BAM files are properly indexed before processing.
    """
    # Generate output file names based on input VCF
    vcf_prefix = os.path.splitext(os.path.basename(vcf_file))[0]
    output_dir = Path(vcf_file).parent
    phased_vcf = os.path.join(output_dir, f"{vcf_prefix}_Phased.vcf")
    phased_vcf_gz = f"{phased_vcf}.gz"
    whatshap_log = os.path.join(output_dir, "whatshap.log")

    # Open log file for WhatsHap output
    with open(whatshap_log, 'w') as log_file:
        # Index the input BAM file
        log_file.write(f"\nIndexing input BAM file: {bam_file}\n")
        index_command = f"samtools index {bam_file}"
        try:
            result = subprocess.run(index_command, shell=True, check=True, capture_output=True, text=True)
            log_file.write("BAM indexing completed successfully\n")
            if result.stderr:
                log_file.write(f"Index stderr: {result.stderr}\n")
        except subprocess.CalledProcessError as e:
            log_file.write(f"Error indexing BAM file: {e.stderr}\n")
            raise

        # Run WhatsHap phase
        phase_command = f"whatshap phase -o {phased_vcf} --reference {reference_fasta} --internal-downsampling 23 --ignore-read-groups {vcf_file} {bam_file}"
        log_file.write("\nRunning WhatsHap phase command:\n")
        log_file.write(f"{phase_command}\n")
        
        try:
            phase_result = subprocess.run(phase_command, shell=True, check=True, capture_output=True, text=True)
            log_file.write("\nWhatsHap Phase stdout:\n")
            log_file.write(phase_result.stdout)
            
            if phase_result.stderr:
                log_file.write("\nWhatsHap Phase stderr:\n")
                log_file.write(phase_result.stderr)
        except subprocess.CalledProcessError as e:
            log_file.write("\nError running WhatsHap phase:\n")
            log_file.write(f"stdout: {e.stdout}\n")
            log_file.write(f"stderr: {e.stderr}\n")
            raise

        # Compress and index the phased VCF
        bgzip_command = f"bgzip -f {phased_vcf}"
        subprocess.run(bgzip_command, shell=True, check=True)
        
        index_command = f"tabix -p vcf {phased_vcf_gz}"
        subprocess.run(index_command, shell=True, check=True)

        # Run WhatsHap haplotag
        haplotag_command = f"whatshap haplotag --tag-supplementary -o {output_bam} --reference {reference_fasta} {phased_vcf_gz} {bam_file} --ignore-read-groups"
        log_file.write("\nRunning WhatsHap haplotag command:\n")
        log_file.write(f"{haplotag_command}\n")
        
        try:
            haplotag_result = subprocess.run(haplotag_command, shell=True, check=True, capture_output=True, text=True)
            log_file.write("\nWhatsHap Haplotag stdout:\n")
            log_file.write(haplotag_result.stdout)
            
            if haplotag_result.stderr:
                log_file.write("\nWhatsHap Haplotag stderr:\n")
                log_file.write(haplotag_result.stderr)
        except subprocess.CalledProcessError as e:
            log_file.write("\nError running WhatsHap haplotag:\n")
            log_file.write(f"stdout: {e.stdout}\n")
            log_file.write(f"stderr: {e.stderr}\n")
            raise

        # Index the output BAM
        log_file.write(f"\nIndexing output BAM file: {output_bam}\n")
        index_bam_command = f"samtools index {output_bam}"
        try:
            subprocess.run(index_bam_command, shell=True, check=True)
            log_file.write("Output BAM indexing completed successfully\n")
        except subprocess.CalledProcessError as e:
            log_file.write(f"Error indexing output BAM file: {e.stderr}\n")
            raise

    return phased_vcf_gz

    
    
def validate_variants(input_bam, variants):
    """
    Validate that the specified variants match the expected patterns in the reads.
    Provides warnings for moderately skewed allele frequencies and errors for severe cases.
    
    Args:
        input_bam: pysam.AlignmentFile object
        variants: list of pysam.VariantRecord objects
    
    Returns:
        tuple: (bool, list, list) - (validation passed, warnings, errors)
    """
    validation_results = []
    warnings = []
    errors = []
    
    for variant in variants:
        ref_count = 0
        alt_count = 0
        unexpected_alleles = Counter()  # Track specific unexpected alleles
        total_covering_reads = 0
        
        # Expected alleles from VCF
        expected_ref = variant.ref.upper()
        expected_alt = variant.alts[0].upper() if variant.alts else ""
        
        # Fetch reads covering this variant position
        for read in input_bam.fetch(variant.chrom, variant.pos - 1, variant.pos):
            allele = get_allele(read, variant)
            if allele == 'ref':
                ref_count += 1
                total_covering_reads += 1
            elif allele == 'alt':
                alt_count += 1
                total_covering_reads += 1
            elif allele:  # Track unexpected alleles
                query_pos = None
                for qpos, rpos in read.get_aligned_pairs():
                    if rpos == variant.pos - 1 and qpos is not None:
                        query_pos = qpos
                        break
                if query_pos is not None:
                    unexpected_base = read.query_sequence[query_pos]
                    unexpected_alleles[unexpected_base] += 1
                    total_covering_reads += 1
        
        if total_covering_reads == 0:
            errors.append(f"No reads cover variant at position {variant.pos}")
            continue
            
        # Calculate percentages
        ref_percentage = (ref_count / total_covering_reads * 100) if total_covering_reads > 0 else 0
        alt_percentage = (alt_count / total_covering_reads * 100) if total_covering_reads > 0 else 0
        unexpected_percentage = sum(unexpected_alleles.values()) / total_covering_reads * 100 if total_covering_reads > 0 else 0
        
        result = {
            'position': variant.pos,
            'total_reads': total_covering_reads,
            'ref_count': ref_count,
            'alt_count': alt_count,
            'ref_percentage': ref_percentage,
            'alt_percentage': alt_percentage,
            'unexpected_alleles': unexpected_alleles,
            'unexpected_percentage': unexpected_percentage,
            'expected_ref': expected_ref,
            'expected_alt': expected_alt
        }
        validation_results.append(result)
        
        # Validation checks
        if total_covering_reads < 30:
            errors.append(
                f"Insufficient coverage at position {variant.pos}:\n"
                f"Total reads: {total_covering_reads} (minimum 30 required)\n"
                "Please verify that the correct variant position was provided."
            )
            continue
        
        # Check for unexpected alleles
        if unexpected_percentage > 0:
            unexpected_details = ", ".join(f"{base}: {count} reads ({count/total_covering_reads*100:.1f}%)" 
                                        for base, count in unexpected_alleles.most_common())
            errors.append(
                f"Unexpected alleles found at position {variant.pos}:\n"
                f"Expected ref: {expected_ref}, Expected alt: {expected_alt}\n"
                f"Unexpected alleles: {unexpected_details}\n"
                "Please verify that the correct variant position and alleles were provided."
            )
        
        # Check allele balance
        if ref_percentage > 95 or alt_percentage > 95:
            # Severe skew - error
            errors.append(
                f"Variant at position {variant.pos} appears to be homozygous:\n"
                f"Reference allele ({expected_ref}): {ref_count} reads ({ref_percentage:.1f}%)\n"
                f"Alternate allele ({expected_alt}): {alt_count} reads ({alt_percentage:.1f}%)\n"
                "This level of skew suggests this position is homozygous. A wrong variant is provided or allele dropout has occured.\n"
                "Please verify that the correct variant position was provided."
            )
        elif ref_percentage > 80 or alt_percentage > 80:
            # Moderate skew - warning
            warnings.append(
                f"Warning: Variant at position {variant.pos} shows skewed allele frequencies:\n"
                f"Reference allele ({expected_ref}): {ref_count} reads ({ref_percentage:.1f}%)\n"
                f"Alternate allele ({expected_alt}): {alt_count} reads ({alt_percentage:.1f}%)\n"
                "While this is acceptable, please verify that the correct variant position was provided."
            )
    
    # Determine if validation passed
    validation_passed = len(errors) == 0
    
    return validation_passed, warnings, errors

def quality_control(input_bam, vcf_file, output_bam):
    """
    Perform quality control on input BAM file and filter reads based on quality metrics.
    """
    if not os.path.exists(input_bam + '.bai'):
        pysam.index(input_bam)
    vcf = pysam.VariantFile(vcf_file)
    variants = list(vcf.fetch())
    if len(variants) != 2:
        raise ValueError("VCF file should contain exactly two variants")

    var1, var2 = variants
    
    input_bam_file = pysam.AlignmentFile(input_bam, "rb")
    
    # Validate variants before proceeding
    variants_valid, warnings, errors = validate_variants(input_bam_file, variants)
    
    print("\n=== Variant Validation ===")
    if variants_valid and not warnings:
        print("Variants were validated successfully")
        print(f"Both variants show expected heterozygous patterns in the reads.")
    elif variants_valid and warnings:
        print("Variants were validated successfully but with warnings:")
        for warning in warnings:
            print(warning)
            print()
    
    if errors:
        print("\nErrors:")
        for error in errors:
            print(error)
            print()
        input_bam_file.close()
        raise ValueError("Variant validation failed. Analysis cannot proceed.")
    
    output_bam_file = pysam.AlignmentFile(output_bam, "wb", template=input_bam_file)

    total_reads = input_bam_file.count()
    spanning_reads = 0
    clean_spanning_reads = 0
    high_quality_reads = 0

    print("\n=== Quality Control ===")

    for read in input_bam_file.fetch():
        if read.reference_start <= var1.pos - 1 and read.reference_end >= var2.pos:
            spanning_reads += 1

            alleles = [get_allele(read, variant) for variant in variants]

            if alleles[0] in ['ref', 'alt'] and alleles[1] in ['ref', 'alt']:
                clean_spanning_reads += 1

                if read.mapping_quality >= 20:
                    high_quality_reads += 1
                    output_bam_file.write(read)

    input_bam_file.close()
    output_bam_file.close()

    print(f"Total reads: {total_reads}")
    print(f"Spanning reads: {spanning_reads} ({spanning_reads/total_reads*100:.2f}%)")
    print(f"Clean spanning reads: {clean_spanning_reads} ({clean_spanning_reads/spanning_reads*100:.2f}% of spanning reads)")
    print(f"High quality spanning reads (MAPQ >= 20): {high_quality_reads} ({high_quality_reads/clean_spanning_reads*100:.2f}% of clean spanning reads)")

    if high_quality_reads > 50:
        print("\n*QC PASSED (>50 high quality spanning reads)*")
    else:
        print("*QC FAILED (<50 high quality spanning reads)*")

    return high_quality_reads

def analyze_reads(bam_file, vcf_file):
    """
    Analyze phased reads to determine variant relationships (Cis/Trans).
    """
    vcf = pysam.VariantFile(vcf_file)
    variants = list(vcf.fetch())
    if len(variants) != 2:
        raise ValueError("VCF file should contain exactly two variants")

    var1, var2 = variants

    bam = pysam.AlignmentFile(bam_file, "rb")

    total_reads = bam.count()
    ref_reads = 0
    alt_reads = 0
    ref_alt_reads = 0
    alt_ref_reads = 0

    for read in bam.fetch():
        alleles = [get_allele(read, variant) for variant in variants]

        if alleles[0] == 'ref' and alleles[1] == 'ref':
            ref_reads += 1
        elif alleles[0] == 'alt' and alleles[1] == 'alt':
            alt_reads += 1
        elif alleles[0] == 'ref' and alleles[1] == 'alt':
            ref_alt_reads += 1
        elif alleles[0] == 'alt' and alleles[1] == 'ref':
            alt_ref_reads += 1

    cis_reads = ref_reads + alt_reads
    trans_reads = ref_alt_reads + alt_ref_reads
    cis_percentage = cis_reads / total_reads * 100 if total_reads > 0 else 0
    trans_percentage = trans_reads / total_reads * 100 if total_reads > 0 else 0
    print("\n=== Result ===")
    print(f"Total high quality spanning reads: {total_reads}")
    print("\nDetailed categorisation of reads:")
    print(f"Reads with ref allele for both variants (Cis): {ref_reads} ({ref_reads/total_reads*100:.2f}%)")
    print(f"Reads with alt allele for both variants (Cis): {alt_reads} ({alt_reads/total_reads*100:.2f}%)")
    print(f"Reads with ref for first, alt for second (Trans): {ref_alt_reads} ({ref_alt_reads/total_reads*100:.2f}%)")
    print(f"Reads with alt for first, ref for second (Trans): {alt_ref_reads} ({alt_ref_reads/total_reads*100:.2f}%)")
    print(f"\nCis reads (both ref or both alt): {cis_reads} ({cis_percentage:.2f}%)")
    print(f"Trans reads (one ref, one alt): {trans_reads} ({trans_percentage:.2f}%)")

    if cis_percentage > trans_percentage:
        print(f"\nCounting the reads determined the phase as Cis \n\nChimeric reads percentage: {trans_percentage:.2f}% ")
    else:
        print(f"\nCounting reads determined the phase as Trans \n\nChimeric reads percentage: {cis_percentage:.2f}% ")

    bam.close()

def parse_vcf_phasing(phased_vcf_file):
    """
    Read the phased VCF file and determine if variants are in Cis or Trans.
    
    Args:
        phased_vcf_file (str): Path to the phased VCF file
    
    Returns:
        str: Interpretation of variant phasing (Cis or Trans)
    """
    vcf = pysam.VariantFile(phased_vcf_file)
    variants = list(vcf.fetch())
    
    if len(variants) != 2:
        raise ValueError("VCF file should contain exactly two variants")
    
    # Extract genotype information
    var1_genotype = variants[0].samples[0]['GT']
    var2_genotype = variants[1].samples[0]['GT']
    
    # Determine phasing
    if var1_genotype == var2_genotype:
        return "Cis"
    else:
        return "Trans"

def main():
    """
    Main function to run the complete analysis pipeline.
    """
    if len(sys.argv) != 3:
        print("Usage: python script.py <input.bam> <input.vcf>")
        sys.exit(1)
    
    bam_file = sys.argv[1]
    vcf_file = sys.argv[2]
    episode = os.path.splitext(os.path.basename(bam_file))[0]
    
    # Create output names
    output_bam = bam_file.removesuffix('.bam') + "_phased.bam"
    output_txt = bam_file.removesuffix('.bam') + "_report.txt"

        # Get BED file path from BAM file path 
    bed_file = bam_file.removesuffix('.bam') + '_coordinate.bed'
    
    # Read amplicon coordinates from BED file
    try:
        with open(bed_file, 'r') as bed:
            line = bed.readline().strip()
            chrom, start, end = line.split()[:3]
            start, end = int(start), int(end)
            amplicon_length = end - start
            region = (chrom, start, end)
    except FileNotFoundError:
        print(f"Error: BED file '{bed_file}' not found")
        sys.exit(1)

    with open(output_txt, 'a') as output_file:
        original_stdout = sys.stdout
        original_stderr = sys.stderr
        sys.stdout = output_file
        sys.stderr = output_file

        try:
            dummy_vcf = vcf_file
            
            # Write amplicon_length and region at the top of the report
            with open(output_txt, 'w') as output_file:  # Use 'w' mode to start fresh
                print(f"Report for: {episode}", file=output_file)
                print(f"{'-'*30}", file=output_file)
                print(f"Amplicon region: {chrom}:{start}-{end}", file=output_file)
                print(f"Amplicon length: {amplicon_length:,} bp", file=output_file)
            
            # Look for the variant calling VCF
            variant_calling_vcf = os.path.join(Path(bam_file).parent._str, f"{episode}.wf_snp.vcf.gz")
            
            if os.path.exists(dummy_vcf) and os.path.exists(variant_calling_vcf):
                write_variant_comparison_results(dummy_vcf, variant_calling_vcf, output_txt)
            else:
                print("Skipping variant comparison due to missing files.")
                print(f"Dummy VCF exists: {os.path.exists(dummy_vcf)}")
                print(f"Variant calling VCF exists: {os.path.exists(variant_calling_vcf)}")
                print(f"Looking for variant calling VCF at: {variant_calling_vcf}")

            clean_span_bam = os.path.join(Path(bam_file).parent._str, f"clean-span-hq.bam")
            high_quality_reads = quality_control(bam_file, vcf_file, clean_span_bam)

            if high_quality_reads > 0:
                phased_vcf_gz = run_whatshap(clean_span_bam, vcf_file, output_bam, REFERENCE_FASTA)

                analyze_reads(output_bam, vcf_file)

                # Unzip WhatsHap phased VCF
                unzip_command = f"gunzip -f {phased_vcf_gz}"
                subprocess.run(unzip_command, shell=True, check=True)

                # Analyze WhatsHap phased VCF
                phased_vcf = os.path.splitext(phased_vcf_gz)[0]
                whatshap_phasing = parse_vcf_phasing(phased_vcf)
                print(" "*40)
                print(f"WhatsHap analysis determined the phase as {whatshap_phasing}")

                # Look for and analyze HapCUT2 phased VCF
                hapcut2_vcf = os.path.join(Path(bam_file).parent._str, f"hap2cut_{episode}.phased.VCF")
                if os.path.exists(hapcut2_vcf):
                    try:
                        hapcut2_phasing = parse_vcf_phasing(hapcut2_vcf)
                        print(f"\nHapCUT2 analysis determined the phase as {hapcut2_phasing}\n")
                    except Exception as e:
                        print(f"Error analyzing HapCUT2 VCF: {e}")
                else:
                    print(f"HapCUT2 phased VCF not found at: {hapcut2_vcf}")
            else:
                print("No high quality spanning reads found. Analysis cannot proceed.")

        except Exception as e:
            print(f"An error occurred: {e}")
            import traceback
            traceback.print_exc()
        finally:
            sys.stdout = original_stdout
            sys.stderr = original_stderr

    print(f"Analysis complete. Results written to {output_txt}")

if __name__ == "__main__":
    main()