#!/usr/bin/env python3

import pysam
import numpy as np
from pathlib import Path
import sys
from collections import defaultdict
from variant_comparison import compare_variants

def calculate_quality_stats(qualities):
    """Calculate mean quality and percentage of bases above Q10"""
    if not qualities:
        return 0, 0
    mean_qual = np.mean(qualities)
    above_q10 = sum(q >= 10 for q in qualities) / len(qualities) * 100
    return mean_qual, above_q10

def calculate_depth_stats(bam, region):
    """Calculate depth statistics for the amplicon region"""
    depths = []
    for pileup_column in bam.pileup(region[0], region[1], region[2], 
                                   truncate=True, min_base_quality=0):
        depths.append(pileup_column.n)
    
    if depths:
        return {
            'median': np.median(depths),
            'mean': np.mean(depths),
            'min': min(depths),
            'max': max(depths)
        }
    return None

def perform_qc_analysis(bam_file, bed_file):
    """
    Perform ONT-specific quality control analysis on BAM file for single amplicon region
    """
    bam_path = Path(bam_file)
    report_file = str(bam_path.parent / f"{bam_path.stem}_report.txt")
    output_bam = str(bam_path.parent / f"{bam_path.stem}_QC.bam")
    
    try:
        # Read the single region from BED file
        with open(bed_file, 'r') as bed:
            line = bed.readline().strip()
            chrom, start, end = line.split()[:3]
            start, end = int(start), int(end)
            amplicon_length = end - start
            region = (chrom, start, end)
        
        # Open input BAM file
        bam = pysam.AlignmentFile(bam_file, "rb")
        out_bam = pysam.AlignmentFile(output_bam, "wb", template=bam)
        
        # Initialize counters and storage
        read_lengths = []
        mapping_quals = []
        base_qualities = []
        read_identities = []
        strand_counts = defaultdict(int)
        passing_reads = []
        
        # Calculate depth statistics
        depth_stats = calculate_depth_stats(bam, region)
        
        # Process reads in the amplicon region
        for read in bam.fetch(region[0], region[1], region[2]):
            # Basic metrics
            read_lengths.append(read.query_length)
            mapping_quals.append(read.mapping_quality)
            
            # Strand information
            strand_counts['forward' if not read.is_reverse else 'reverse'] += 1
            
            # Base qualities
            if read.query_qualities is not None:
                mean_qual, above_q10 = calculate_quality_stats(read.query_qualities)
                base_qualities.append(mean_qual)
            
            # Calculate read identity
            aligned_pairs = read.get_aligned_pairs(with_seq=True)
            matches = sum(1 for _, _, base in aligned_pairs if base and base.isupper())
            identity = matches / len(aligned_pairs) if aligned_pairs else 0
            read_identities.append(identity)
            
            # Store reads that pass QC criteria
            if (read.query_length >= 1000 and 
                read.mapping_quality >= 20 and 
                identity >= 0.8 and
                (not read.query_qualities or 
                 calculate_quality_stats(read.query_qualities)[0] >= 7)):
                passing_reads.append(read)
        
        total_reads = len(read_lengths)
        passing_qc_reads = len(passing_reads)
        min_depth = depth_stats['min'] if depth_stats else 0
        
        # Write QC report
        with open(report_file, 'w') as report:
            report.write(f"Report for {bam_path.stem}\n")
            report.write("-" * 30 + "\n")
            
            # Add amplicon region, length, and passing QC reads to the top
            report.write(f"Amplicon region: {chrom}:{start}-{end}\n")
            report.write(f"Amplicon length: {amplicon_length:,} bp\n")
            report.write(f"Total reads: {total_reads}\n")
            report.write(f"Passing QC reads: {passing_qc_reads} ({passing_qc_reads/total_reads*100:.2f}%)\n")
            
            # Add QC check messages
            if passing_qc_reads > 50 and min_depth > 50:
                report.write("\n*QC PASSED (>50 HQ reads and minimum depth >50x)*\n")
            elif passing_qc_reads > 50 and min_depth < 50:
                report.write("\n*QC FAILED (Minimum depth <50x)*\n")
            elif passing_qc_reads < 50 and min_depth > 50:
                report.write("\n*QC FAILED (HQ reads <50)*\n")
            else:
                report.write("\n*QC FAILED (<50 HQ reads and minimum depth <50x)*\n")
            
            report.write("\n")
            
            # Check for variant calling VCF
            episode = bam_path.stem
            variant_calling_vcf = bam_path.parent / f"{episode}.wf_snp.vcf.gz"
            dummy_vcf = bam_path.parent / f"{episode}.vcf"
            
            # Add variant comparison results if files exist
            if dummy_vcf.exists() and variant_calling_vcf.exists():
                report.write("\n=== Variant Calling ===\n")
                report.write("User provided variant:\n")
                # Read dummy VCF variants and show only essential columns
                with open(str(dummy_vcf), 'r') as dummy_f:
                    for line in dummy_f:
                        if not line.startswith('#'):
                            parts = line.strip().split('\t')
                            if len(parts) >= 5:
                                chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
                                report.write(f"{chrom}\t{pos}\t{ref}\t{alt}\n")
                
                # Create a set of dummy variants for matching
                dummy_variants_set = set()
                with open(str(dummy_vcf), 'r') as dummy_f:
                    for line in dummy_f:
                        if not line.startswith('#'):
                            parts = line.strip().split('\t')
                            if len(parts) >= 5:
                                chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
                                dummy_variants_set.add((chrom, int(pos), ref, alt))
                
                report.write("\nVariants called from the amplicon by Clair3:\n")
                report.write("CHROM\tPOS\tREF\tALT\tQUAL\tGT\n")
                # Read variant calling VCF variants
                vcf = pysam.VariantFile(str(variant_calling_vcf))
                for record in vcf.fetch():
                    if record.alts and record.alts[0] and record.alts[0] != ".":
                        qual = f"{record.qual:.1f}" if record.qual is not None else "."
                        gt = record.samples[0]['GT']
                        gt_str = "/".join(str(x) for x in gt) if gt is not None else "."
                        
                        # Check if this variant matches any in the dummy VCF
                        is_matched = (record.chrom, record.pos, record.ref, record.alts[0]) in dummy_variants_set
                        
                        # Add "<-" for matched variants
                        indicator = "<-" if is_matched else " "
                        report.write(f"{record.chrom}\t{record.pos}\t{record.ref}\t{record.alts[0]}\t{qual}\t{gt_str} {indicator}\n")
                
                # Add variant matching results
                matched_variants = compare_variants(str(dummy_vcf), str(variant_calling_vcf))
                report.write("\nVariant Matching:\n")
                for i, variant in enumerate(matched_variants, 1):
                    if variant:
                        report.write(f"Variant {i}: {variant[0]} {variant[1]} {variant[2]} > {variant[3]} - FOUND\n")
                    else:
                        report.write(f"Variant {i}: NOT FOUND in variant calling VCF\n")
            
            report.write(f"\n==Quality control details==\n")
            
            if depth_stats:
               report.write(f"Depth Statistics:\n")
               report.write(f"Median depth: {depth_stats['median']:.0f}X\n")
               report.write(f"Minimum depth: {depth_stats['min']:.0f}X\n")
               report.write(f"Maximum depth: {depth_stats['max']:.0f}X\n")
            
            report.write(f"\nRead Length Distribution:\n")
            report.write(f"Median length: {np.median(read_lengths):.0f}\n")
            report.write(f"Mean length: {np.mean(read_lengths):.0f}\n")
            report.write(f"Length N50: {calculate_n50(read_lengths):.0f}\n\n")
            
            report.write(f"Mapping Quality:\n")
            report.write(f"Mean MAPQ: {np.mean(mapping_quals):.0f}\n")
            report.write(f"Reads with MAPQ≥20: {sum(q >= 20 for q in mapping_quals)} ")
            report.write(f"({sum(q >= 20 for q in mapping_quals)/total_reads*100:.0f}%)\n\n")
            
            if base_qualities:
                report.write(f"Base Quality:\n")
                report.write(f"Mean base quality: {np.mean(base_qualities):.0f}\n")
                report.write(f"Median base quality: {np.median(base_qualities):.0f}\n\n")
            
            report.write(f"Alignment Statistics:\n")
            report.write(f"Mean identity: {np.mean(read_identities)*100:.0f}%\n")
            report.write(f"Reads ≥80% identity: {sum(i >= 0.8 for i in read_identities)} ")
            report.write(f"({sum(i >= 0.8 for i in read_identities)/total_reads*100:.0f}%)\n\n")
            
            report.write(f"Strand Distribution:\n")
            report.write(f"Forward: {strand_counts['forward']} ({strand_counts['forward']/total_reads*100:.0f}%)\n")
            report.write(f"Reverse: {strand_counts['reverse']} ({strand_counts['reverse']/total_reads*100:.0f}%)\n")
            
        # Write passing reads to output BAM file
        for read in passing_reads:
            out_bam.write(read)
            
        print(f"Report has been saved to: {report_file}")
        
        # Close and index BAM file
        out_bam.close()
        bam.close()
        pysam.index(output_bam)
            
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)
    finally:
        if 'bam' in locals():
            bam.close()
        if 'out_bam' in locals():
            out_bam.close()

def calculate_n50(lengths):
    """Calculate N50 for a list of lengths"""
    if not lengths:
        return 0
    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)
    running_sum = 0
    for length in sorted_lengths:
        running_sum += length
        if running_sum >= total_length / 2:
            return length
    return 0

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <bam_file> <bed_file>")
        sys.exit(1)
    
    bam_file = sys.argv[1]
    bed_file = sys.argv[2]
    
    if not Path(bam_file).exists():
        print(f"Error: BAM file '{bam_file}' not found")
        sys.exit(1)
    if not Path(bed_file).exists():
        print(f"Error: BED file '{bed_file}' not found")
        sys.exit(1)
        
    perform_qc_analysis(bam_file, bed_file)