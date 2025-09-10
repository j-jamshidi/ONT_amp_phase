import pysam
import os
import sys
from typing import List, Tuple

def compare_variants(dummy_vcf: str, variant_calling_vcf: str) -> List[Tuple[str, int, str, str]]:
    """
    Compare variants between dummy VCF and variant calling VCF.
    
    Args:
        dummy_vcf (str): Path to the dummy VCF file
        variant_calling_vcf (str): Path to the variant calling VCF file
    
    Returns:
        List of variants found in both files, with details about matching status
    """
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
    """
    Write variant comparison results to the report file.
    
    Args:
        dummy_vcf (str): Path to the dummy VCF file
        variant_calling_vcf (str): Path to the variant calling VCF file
        output_file (str): Path to the output report file
    """
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
        
        with open(output_file, 'a') as f:
            # Create a set of dummy variants for matching
            dummy_variants_set = set()
            with open(dummy_vcf, 'r') as dummy_f:
                for line in dummy_f:
                    if not line.startswith('#'):
                        parts = line.strip().split('\t')
                        if len(parts) >= 5:
                            chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
                            dummy_variants_set.add((chrom, int(pos), ref, alt))

            f.write("===Variant Calling===\n")           
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
                    
                    # Add "<-" for matched variants
                    indicator = "<-" if is_matched else " "
                    f.write(f"{record.chrom}\t{record.pos}\t{record.ref}\t{record.alts[0]}\t{qual}\t{gt_str} {indicator} \n")
            
            f.write("\nVariant Matching Results:\n")
            for i, variant in enumerate(matched_variants, 1):
                if variant:
                    f.write(f"Variant {i}: {variant[0]} {variant[1]} {variant[2]} > {variant[3]} - FOUND\n")
                else:
                    f.write(f"Variant {i}: NOT FOUND in variant calling VCF\n")
    
    except Exception as e:
        print(f"Error during variant comparison: {e}")
        import traceback
        traceback.print_exc()

def main():
    """
    Main function to run variant comparison standalone.
    """
    if len(sys.argv) != 4:
        print("Usage: python variant_comparison.py <dummy_vcf> <variant_calling_vcf> <output_file>")
        sys.exit(1)
    
    dummy_vcf = sys.argv[1]
    variant_calling_vcf = sys.argv[2]
    output_file = sys.argv[3]
    
    write_variant_comparison_results(dummy_vcf, variant_calling_vcf, output_file)
    print(f"Variant comparison completed. Results written to {output_file}")

if __name__ == "__main__":
    main()