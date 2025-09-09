import pysam
import argparse
import sys
import os
from collections import Counter
from typing import Optional, Tuple, List
from pathlib import Path

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

def main():
    """
    Main function to run quality control and create clean spanning BAM file.
    """
    if len(sys.argv) != 3:
        print("Usage: python phasing_variants_qc.py <input.bam> <input.vcf>")
        sys.exit(1)
    
    bam_file = sys.argv[1]
    vcf_file = sys.argv[2]
    episode = os.path.splitext(os.path.basename(bam_file))[0]
    output_dir = os.path.dirname(bam_file)
    
    # Create output files
    clean_span_bam = os.path.join(output_dir, "clean-span-hq.bam")
    report_file = os.path.join(output_dir, f"{episode}_report.txt")
    
    # Get BED file path and read amplicon coordinates
    base_name = bam_file[:-4] if bam_file.endswith('.bam') else bam_file
    bed_file = base_name + '_coordinate.bed'
    try:
        with open(bed_file, 'r') as bed:
            line = bed.readline().strip()
            chrom, start, end = line.split()[:3]
            start, end = int(start), int(end)
            amplicon_length = end - start
    except FileNotFoundError:
        print(f"Error: BED file '{bed_file}' not found")
        sys.exit(1)
    
    # Write initial report info
    with open(report_file, 'w') as f:
        print(f"Report for: {episode}", file=f)
        print(f"{'-'*30}", file=f)
        print(f"Amplicon region: {chrom}:{start}-{end}", file=f)
        print(f"Amplicon length: {amplicon_length:,} bp", file=f)
    
    # Append QC results to report
    with open(report_file, 'a') as f:
        original_stdout = sys.stdout
        sys.stdout = f
        
        try:
            high_quality_reads = quality_control(bam_file, vcf_file, clean_span_bam)
            
            # Index the output BAM
            pysam.index(clean_span_bam)
            
        except Exception as e:
            print(f"Error during quality control: {e}")
            sys.exit(1)
        finally:
            sys.stdout = original_stdout
    
    print(f"Quality control completed. Results written to {report_file}")
    print(f"Clean spanning BAM created: {clean_span_bam}")

if __name__ == "__main__":
    main()