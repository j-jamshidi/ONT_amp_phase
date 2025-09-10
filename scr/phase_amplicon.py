import pysam
import os
import sys
from pathlib import Path
from phasing_variants_qc import get_allele


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
        phased_vcf_file (str): Path to the phased VCF file (can be gzipped)
    
    Returns:
        str: Interpretation of variant phasing (Cis or Trans)
    """
    vcf = pysam.VariantFile(phased_vcf_file)  # pysam automatically handles gzipped files
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

def analyze_phasing_results(episode, base_dir):
    """
    Analyze all available phasing VCF files and report results.
    
    Args:
        episode (str): Sample episode name
        base_dir (str): Base directory containing VCF files
    """
    # Analyze WhatsHap phased VCF (gzipped)
    phased_vcf_gz = os.path.join(base_dir, f"{episode}_Phased.vcf.gz")
    if os.path.exists(phased_vcf_gz):
        whatshap_phasing = parse_vcf_phasing(phased_vcf_gz)
        print(" "*40)
        print(f"WhatsHap analysis determined the phase as {whatshap_phasing}")
    else:
        print("WhatsHap phased VCF not found.")

    # Look for and analyze HapCUT2 phased VCF
    hapcut2_vcf = os.path.join(base_dir, f"hap2cut_{episode}.phased.VCF")
    if os.path.exists(hapcut2_vcf):
        try:
            hapcut2_phasing = parse_vcf_phasing(hapcut2_vcf)
            print(f"\nHapCUT2 analysis determined the phase as {hapcut2_phasing}\n")
        except Exception as e:
            print(f"Error analyzing HapCUT2 VCF: {e}")
    else:
        print(f"HapCUT2 phased VCF not found at: {hapcut2_vcf}")

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
    base_name = bam_file[:-4] if bam_file.endswith('.bam') else bam_file
    output_bam = base_name + "_phased.bam"
    output_txt = base_name + "_report.txt"

    with open(output_txt, 'a') as output_file:
        original_stdout = sys.stdout
        original_stderr = sys.stderr
        sys.stdout = output_file
        sys.stderr = output_file

        try:
            dummy_vcf = vcf_file
            


            # Quality control is now handled by separate script
            clean_span_bam = os.path.join(os.path.dirname(bam_file), "clean-span-hq.bam")
            
            if os.path.exists(clean_span_bam):
                # WhatsHap is now run from the shell script
                analyze_reads(output_bam, vcf_file)
                
                # Analyze all phasing results
                analyze_phasing_results(episode, os.path.dirname(bam_file))
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