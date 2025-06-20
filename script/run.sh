#!/bin/bash

# Script: run.sh
# Description: Pipeline for processing and phasing ONT amplicon data.
# Works for phasing when two variants are provided, or localisation/QC when one or none is provided.
#===============================================================================
# Configuration and Setup
#===============================================================================
RUNID=$1

#Change "/path/to/directory" to your working directory
BASEDIR="/path/to/directory/ONT/Runs/${RUNID}"

WORKDIR="${BASEDIR}/result"

# Reference and tool paths. Change as required
REFERENCE="/path/to/reference/directory/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
#CLAIR3_PATH="/EBSDataDrive/ONT/Clair3" OR use docker (in the current setting)
HAPCUT2_PATH="/path/to/HapCUT2/build"
SCRIPT_PATH="/path/to/directory/ONT/script"
#===============================================================================
# Functions
#===============================================================================

# Logging function
log() {
    echo -e "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Prepare VCF file
prepare_vcf() {
    local episode=$1
    local barcode=$2
    local variant1=$3
    local variant2=$4

    cp "${SCRIPT_PATH}/dummy.vcf" "${WORKDIR}/${barcode}/${episode}.vcf"
    sed -i -e "s/sample/${episode}/g" "${WORKDIR}/${barcode}/${episode}.vcf"

    # Process first variant
    variant1=$(echo "${variant1}" | tr ':>\t ' '-' | sed 's/--*/-/g')
    IFS='-' read -r -a array1 <<< "${variant1}"
    sed -i -e "s/chrA/${array1[0]}/g" "${WORKDIR}/${barcode}/${episode}.vcf"
    sed -i -e "s/POS1/${array1[1]}/g" "${WORKDIR}/${barcode}/${episode}.vcf"
    sed -i -e "s/REF1/${array1[2]}/g" "${WORKDIR}/${barcode}/${episode}.vcf"
    sed -i -e "s/ALT1/${array1[3]}/g" "${WORKDIR}/${barcode}/${episode}.vcf"
    unset array1

    if [[ "$variant2" =~ chr ]]; then
        # Process second variant
        variant2=$(echo "${variant2}" | tr ':>\t ' '-' | sed 's/--*/-/g')
        IFS='-' read -r -a array2 <<< "${variant2}"
        sed -i -e "s/chrB/${array2[0]}/g" "${WORKDIR}/${barcode}/${episode}.vcf"
        sed -i -e "s/POS2/${array2[1]}/g" "${WORKDIR}/${barcode}/${episode}.vcf"
        sed -i -e "s/REF2/${array2[2]}/g" "${WORKDIR}/${barcode}/${episode}.vcf"
        sed -i -e "s/ALT2/${array2[3]}/g" "${WORKDIR}/${barcode}/${episode}.vcf"
        unset array2

        # Sort VCF file
        grep -v "^#" "${WORKDIR}/${barcode}/${episode}.vcf" | sort -k1,1 -k2,2n > "${WORKDIR}/${barcode}/${episode}_sorted.vcf"
        grep "^#" "${WORKDIR}/${barcode}/${episode}.vcf" > "${WORKDIR}/${barcode}/${episode}_header.vcf"
        cat "${WORKDIR}/${barcode}/${episode}_header.vcf" "${WORKDIR}/${barcode}/${episode}_sorted.vcf" > "${WORKDIR}/${barcode}/${episode}.vcf"
        rm "${WORKDIR}/${barcode}/${episode}_header.vcf" "${WORKDIR}/${barcode}/${episode}_sorted.vcf"
    else
        # Remove second variant line for single variant case
        sed -i '/chrB/d' "${WORKDIR}/${barcode}/${episode}.vcf"
    fi
}

# Merge BAM files
merge_bam_files() {
    local barcode=$1
    local episode=$2
    local coordinate=$3

    log "Merging BAM files for ${barcode}, sample ${episode}..."
         samtools merge -@ 6 -u - "${BASEDIR}/bam_pass/${barcode}"/*bam 2>/dev/null | \
         samtools sort -@ 6 -o "${WORKDIR}/${barcode}/temp.bam" 2>/dev/null && \
         samtools index "${WORKDIR}/${barcode}/temp.bam" 2>/dev/null && \
         samtools view -@ 6 -b "${WORKDIR}/${barcode}/temp.bam" "${coordinate}" > "${WORKDIR}/${barcode}/${episode}.bam" 2>/dev/null && \
         samtools index "${WORKDIR}/${barcode}/${episode}.bam" 2>/dev/null && \
     rm "${WORKDIR}/${barcode}/temp"*
log "BAM files merged!"
}

# Run Clair3 for variant calling
run_clair3() {
    local barcode=$1
    local episode=$2

    log "Running variant calling with Clair3 for ${barcode}, sample ${episode}..."
    docker run --rm --platform linux/amd64 \
    -v "${WORKDIR}:/data" \
    -v "$(dirname "${REFERENCE}"):/refs" \
    hkubal/clair3:latest \
    /opt/bin/run_clair3.sh \
    --bam_fn="/data/${barcode}/${episode}.bam" \
    --bed_fn="/data/${barcode}/${episode}_coordinate.bed" \
    --ref_fn="/refs/$(basename "${REFERENCE}")" \
    --threads=6 \
    --platform=ont \
    --model_path="/opt/models/r1041_e82_400bps_sup_v500" \
    --sample_name="${episode}" \
    --use_whatshap_for_final_output_phasing \
    --enable_phasing \
    --remove_intermediate_dir \
    --var_pct_full=1 \
    --ref_pct_full=1 \
    --var_pct_phasing=1 \
    --output="/data/${barcode}/variant_calling_output" \
    --min_coverage=20 >/dev/null 2>&1

    log "Clair3 analysis finished!"
}

# Run HapCUT2 for phasing
run_hapcut2() {
    local barcode=$1
    local episode=$2

    log "Running HapCUT2 for ${barcode}, ${episode}..."
    "${HAPCUT2_PATH}/extractHAIRS" \
        --ont 1 \
        --bam "${WORKDIR}/${barcode}/${episode}.bam" \
        --VCF "${WORKDIR}/${barcode}/${episode}.vcf" \
        --out "${WORKDIR}/${barcode}/fragment_${episode}" \
        --indels 1 \
        --ref "${REFERENCE}" > "${WORKDIR}/${barcode}/HapCUT2.log" 2>&1

    "${HAPCUT2_PATH}/HAPCUT2" \
        --fragments "${WORKDIR}/${barcode}/fragment_${episode}" \
        --VCF "${WORKDIR}/${barcode}/${episode}.vcf" \
        --output "${WORKDIR}/${barcode}/hap2cut_${episode}" >> "${WORKDIR}/${barcode}/HapCUT2.log" 2>&1
    log "HapCUT2 analysis finished!"
}

# Prepare input file
prepare_input_file() {
    # Remove header from CSV file and create info file
    tail -n +2 "${BASEDIR}/sample_sheet.csv" > "${BASEDIR}/temp.info"

    # Process the Barcode column to add leading zeros and 'barcode' prefix
    # Convert Episode and EpisodeWES to uppercase
    awk -F',' 'BEGIN {OFS=","} {
        gsub(/[[:space:]]*$/, "", $0)
        if ($2 ~ /^[0-9]+$/) {
            $2 = sprintf("barcode%02d", $2)
        }
        $3 = toupper($3)
        print $0
    }' "${BASEDIR}/temp.info" > "${BASEDIR}/${RUNID}.info"

    # Ensure there's exactly one newline at the end of the file
    sed -i -e '$a\' "${BASEDIR}/${RUNID}.info"

    # Clean up temporary file
    rm "${BASEDIR}/temp.info"
}

# Main processing loop
process_samples() {
    while IFS=, read -r Batch Barcode Episode Coordinate Variant1 Variant2 remainder; do
            log "Processing ${Barcode}..."
            mkdir -p "${WORKDIR}/${Barcode}"

            # Clean up Coordinate by removing trailing spaces
            Coordinate=$(echo "${Coordinate}" | sed 's/[[:space:]]*$//')

            # VCF File Preparation
            if [[ ! "$Variant1" =~ chr ]]; then
                log "No variant is provided, continuing with basecalling and QC..."
            elif [[ ! "$Variant2" =~ chr ]]; then
                log "Only one variant is provided, continuing with basecalling and QC..."
                prepare_vcf "$Episode" "$Barcode" "$Variant1" ""
            else
                log "Two variants are provided, continuing with basecalling, QC, and phasing..."
                prepare_vcf "$Episode" "$Barcode" "$Variant1" "$Variant2"
            fi

            # BAM File Processing
            merge_bam_files "$Barcode" "$Episode" "$Coordinate"

            # Making coordinate bed file
            echo -e "${Coordinate}" | tr -d ' ' | tr ':' '\t' | sed 's/-/\t/g' > "${WORKDIR}/${Barcode}/${Episode}_coordinate.bed"

            # Variant Calling with Clair3
            current_dir=$PWD
            cd "${WORKDIR}/${Barcode}" || exit
            run_clair3 "$Barcode" "$Episode"

            # Copy output files
            cp "${WORKDIR}/${Barcode}/variant_calling_output/merge_output.vcf.gz" "${Episode}.wf_snp.vcf.gz"
            cp "${WORKDIR}/${Barcode}/variant_calling_output/merge_output.vcf.gz.tbi" "${Episode}.wf_snp.vcf.gz.tbi"

            # HapCUT2 Phasing (only for two variants)
            if [[ "$Variant1" =~ chr ]] && [[ "$Variant2" =~ chr ]]; then
                cd "$current_dir" || exit
                run_hapcut2 "$Barcode" "$Episode"
            fi

            # Final Analysis and Cleanup
            log "Analyzing the reads and writing the results for ${Barcode}, ${Episode}..."
            if [[ ! "$Variant1" =~ chr ]] || [[ ! "$Variant2" =~ chr ]]; then
                python "${SCRIPT_PATH}/basecalling_QC_amplicon.py" \
                    "${WORKDIR}/${Barcode}/${Episode}.bam" \
                    "${WORKDIR}/${Barcode}/${Episode}_coordinate.bed"
            else
                python "${SCRIPT_PATH}/basecalling_phasing_amplicon.py" \
                    "${WORKDIR}/${Barcode}/${Episode}.bam" \
                    "${WORKDIR}/${Barcode}/${Episode}.vcf"
            fi

            # Cleanup temporary files
            rm -rf "${WORKDIR}/${Barcode}/work"
            if [[ "$Variant1" =~ chr ]] && [[ "$Variant2" =~ chr ]]; then
                rm "${WORKDIR}/${Barcode}/${Episode}.vcf"
                rm "${WORKDIR}/${Barcode}/hap2cut_${Episode}"
                rm "${WORKDIR}/${Barcode}/fragment_${Episode}"
                rm "${WORKDIR}/${Barcode}/clean-span-hq.bam"*
            fi

    done < "${BASEDIR}/${RUNID}.info"

     log "Done!\n"  
}

#===============================================================================
# Main Script Execution
#===============================================================================

# Prepare input file
prepare_input_file

# Process samples
process_samples

log "Pipeline execution completed!"
