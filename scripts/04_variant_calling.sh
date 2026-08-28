#!/usr/bin/env bash

###############################################################################
## QUINOA VARIANT CALLING AND SNP FILTRATION PIPELINE
###############################################################################
## Variant callingworkflow adapted after Patiranage et al. (2022).
## Reference: Chenopodium quinoa QQ74 version 2 (CoGe ID 60716).
## Cohort: 91 accessions; joint genotyping performed for 18 chromosomes.
##
## This file retains the final production logic from the original working notes.
## Exploratory one-sample commands, failed alternatives, installation commands,
## interactive SSH commands, and data-transfer commands have been removed.
##
## Run one stage at a time:
##   bash organized_CAU_quinoa_variant_calling_pipeline.sh <stage>
##
## Stages:
##   download | trim | reference | align | read-groups | index-bam | call-gvcf
##   compress-gvcf | sample-map | joint-genotype | extract-snps
##   remove-missing-ranks | hard-filter | select-pass | missingness | merge
##
###############################################################################

set -euo pipefail
shopt -s nullglob

###############################################################################
## 0. SOFTWARE AND PROJECT PATHS
###############################################################################

JAVA_HOME="${JAVA_HOME:-/variant_calling/Tools/jdk-23.0.1}"
SRA_BIN="${SRA_BIN:-/variant_calling/Tools/sratoolkit.3.1.11-ubuntu64/bin}"
TRIMMOMATIC_HOME="${TRIMMOMATIC_HOME:-/variant_calling/Tools/Trimmomatic-0.39}"
BWA="${BWA:-//variant_calling/Tools/bwa/bwa}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-/variant_calling/Tools/samtools-1.18}"
GATK_JAR="${GATK_JAR:-/variant_calling/Tools/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar}"

export JAVA_HOME
export PATH="${JAVA_HOME}/bin:${SRA_BIN}:${SAMTOOLS_BIN}:${PATH}"

PROJECT_DIR="${PROJECT_DIR:-/vc}"
REF="${REF:-/Chenopodium_quinoa.fasta}"
RAW_DIR="${RAW_DIR:-/variant_calling/raw_reads_sra/raw_reads}"
TRIM_DIR="${TRIM_DIR:-${PROJECT_DIR}/trimmed_reads}"
SORTED_BAM_DIR="${SORTED_BAM_DIR:-${PROJECT_DIR}/sorted_bam}"
DEDUP_BAM_DIR="${DEDUP_BAM_DIR:-${PROJECT_DIR}/dup_bam}"
GVCF_DIR="${GVCF_DIR:-${PROJECT_DIR}/ind_vcf}"
JOINT_DIR="${JOINT_DIR:-${PROJECT_DIR}/chr_vcf}"
SNP_DIR="${SNP_DIR:-${JOINT_DIR}/snps}"
FINAL_CHR_DIR="${FINAL_CHR_DIR:-${SNP_DIR}/passed_snps}"

SRA_LIST="${SRA_LIST:-${PROJECT_DIR}/SRA_ids.txt}"
SAMPLE_LIST="${SAMPLE_LIST:-${PROJECT_DIR}/samples_list.txt}"
CHROMOSOME_MAP="${CHROMOSOME_MAP:-${PROJECT_DIR}/chromosome_intervals.tsv}"
SAMPLE_NAME_MAP="${SAMPLE_NAME_MAP:-${GVCF_DIR}/sample_name_map.txt}"

THREADS="${THREADS:-16}"
DB_THREADS="${DB_THREADS:-32}"
HETEROZYGOSITY="${HETEROZYGOSITY:-0.005}"
JAVA_HEAP="${JAVA_HEAP:-100g}"

gatk() {
    java "-Xmx${JAVA_HEAP}" -jar "${GATK_JAR}" "$@"
}

require_file() {
    [[ -s "$1" ]] || { echo "ERROR: required file is absent or empty: $1" >&2; exit 1; }
}

sample_count() {
    awk 'NF && $1 !~ /^#/' "${SAMPLE_LIST}" | wc -l
}

chromosome_count() {
    awk 'NF >= 2 && $1 !~ /^#/' "${CHROMOSOME_MAP}" | wc -l
}

###############################################################################
## 1. DOWNLOAD PAIRED-END READS FROM NCBI SRA
###############################################################################

download_reads() {
    require_file "${SRA_LIST}"
    mkdir -p "${RAW_DIR}"

    while IFS= read -r accession; do
        accession="${accession%$'\r'}"
        [[ -n "${accession}" && "${accession}" != \#* ]] || continue
        fastq-dump --split-files --gzip --outdir "${RAW_DIR}" "${accession}"
    done < "${SRA_LIST}"
}

###############################################################################
## 2. TRIM READS WITH TRIMMOMATIC 0.39
###############################################################################
## Final parameters: LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:50.
## Only paired outputs are used for alignment; unpaired outputs are retained.

trim_reads() {
    mkdir -p "${TRIM_DIR}"

    local r1 r2 sample
    for r1 in "${RAW_DIR}"/*_1.fastq.gz; do
        sample="$(basename "${r1}" _1.fastq.gz)"
        r2="${RAW_DIR}/${sample}_2.fastq.gz"
        [[ -f "${r2}" ]] || { echo "WARNING: missing R2 for ${sample}; skipped" >&2; continue; }

        java -jar "${TRIMMOMATIC_HOME}/trimmomatic-0.39.jar" PE \
            -threads "${THREADS}" \
            "${r1}" "${r2}" \
            "${TRIM_DIR}/${sample}_1.paired.fastq.gz" \
            "${TRIM_DIR}/${sample}_1.unpaired.fastq.gz" \
            "${TRIM_DIR}/${sample}_2.paired.fastq.gz" \
            "${TRIM_DIR}/${sample}_2.unpaired.fastq.gz" \
            LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:50
    done
}

###############################################################################
## 3. PREPARE THE QQ74 VERSION 2 REFERENCE GENOME
###############################################################################

prepare_reference() {
    require_file "${REF}"
    "${BWA}" index "${REF}"
    samtools faidx "${REF}"
    gatk CreateSequenceDictionary -R "${REF}"
}

###############################################################################
## 4. ALIGN READS WITH BWA 0.7.18; FILTER AND SORT WITH SAMTOOLS 1.18
###############################################################################
## The original executed commands additionally used BWA-MEM -M -k 30 and
## retained alignments with mapping quality >=30 via samtools view -q 30.

align_reads() {
    require_file "${REF}"
    mkdir -p "${SORTED_BAM_DIR}"

    local r1 r2 sample
    for r1 in "${TRIM_DIR}"/*_1.paired.fastq.gz; do
        sample="$(basename "${r1}" _1.paired.fastq.gz)"
        r2="${TRIM_DIR}/${sample}_2.paired.fastq.gz"
        [[ -f "${r2}" ]] || { echo "WARNING: missing paired R2 for ${sample}; skipped" >&2; continue; }

        "${BWA}" mem -M -k 30 -t "${THREADS}" "${REF}" "${r1}" "${r2}" |
            samtools view -@ "${THREADS}" -b -h -q 30 - |
            samtools sort -@ "${THREADS}" -o "${SORTED_BAM_DIR}/${sample}.sorted.bam" -
    done
}

###############################################################################
## 5. ADD READ GROUPS AND MARK DUPLICATES WITH GATK/PICARD 4.6.1.0
###############################################################################

add_read_groups_and_mark_duplicates() {
    mkdir -p "${DEDUP_BAM_DIR}"

    local bam sample rg_bam
    for bam in "${SORTED_BAM_DIR}"/*.sorted.bam; do
        sample="$(basename "${bam}" .sorted.bam)"
        rg_bam="${DEDUP_BAM_DIR}/${sample}.rgroup.bam"

        gatk AddOrReplaceReadGroups \
            --INPUT "${bam}" \
            --OUTPUT "${rg_bam}" \
            --SORT_ORDER coordinate \
            --RGSM "${sample}" \
            --RGPU none \
            --RGID 1 \
            --RGLB lib1 \
            --RGPL Illumina

        gatk MarkDuplicates \
            --INPUT "${rg_bam}" \
            --OUTPUT "${DEDUP_BAM_DIR}/${sample}.rmdup.bam" \
            --METRICS_FILE "${DEDUP_BAM_DIR}/${sample}.metrics.txt"
    done
}

###############################################################################
## 6. INDEX DEDUPLICATED BAM FILES WITH SAMTOOLS 1.18
###############################################################################

index_bams() {
    local bam
    for bam in "${DEDUP_BAM_DIR}"/*.rmdup.bam; do
        samtools index -@ "${THREADS}" "${bam}"
    done
}

###############################################################################
## 7. CALL PER-SAMPLE GVCFs WITH GATK HAPLOTYPECALLER 4.6.1.0
###############################################################################
## Final recorded setting: --emit-ref-confidence GVCF --heterozygosity 0.005.

call_gvcfs() {
    require_file "${SAMPLE_LIST}"
    [[ "$(sample_count)" -eq 91 ]] || {
        echo "ERROR: expected 91 samples in ${SAMPLE_LIST}; found $(sample_count)" >&2
        exit 1
    }
    mkdir -p "${GVCF_DIR}"

    local sample bam
    while IFS= read -r sample; do
        sample="${sample%$'\r'}"
        [[ -n "${sample}" && "${sample}" != \#* ]] || continue
        bam="${DEDUP_BAM_DIR}/${sample}.rmdup.bam"
        require_file "${bam}"
        require_file "${bam}.bai"

        gatk HaplotypeCaller \
            -R "${REF}" \
            -I "${bam}" \
            --emit-ref-confidence GVCF \
            --heterozygosity "${HETEROZYGOSITY}" \
            -O "${GVCF_DIR}/${sample}.snps.indels.g.vcf"
    done < "${SAMPLE_LIST}"
}

###############################################################################
## 8. BGZIP-COMPRESS AND TABIX-INDEX GVCFs
###############################################################################

compress_gvcfs() {
    local vcf
    for vcf in "${GVCF_DIR}"/*.snps.indels.g.vcf; do
        bgzip -@ "${THREADS}" -c "${vcf}" > "${vcf}.gz"
        tabix -f -p vcf "${vcf}.gz"
    done
}

###############################################################################
## 9. CREATE THE GENOMICSDBIMPORT SAMPLE MAP
###############################################################################
## GATK sample maps contain two columns: sample name and absolute gVCF path.

create_sample_map() {
    require_file "${SAMPLE_LIST}"
    : > "${SAMPLE_NAME_MAP}"

    local sample gvcf
    while IFS= read -r sample; do
        sample="${sample%$'\r'}"
        [[ -n "${sample}" && "${sample}" != \#* ]] || continue
        gvcf="${GVCF_DIR}/${sample}.snps.indels.g.vcf.gz"
        require_file "${gvcf}"
        require_file "${gvcf}.tbi"
        printf '%s\t%s\n' "${sample}" "${gvcf}" >> "${SAMPLE_NAME_MAP}"
    done < "${SAMPLE_LIST}"
}

###############################################################################
## 10. JOINT GENOTYPING BY CHROMOSOME
###############################################################################
## GenomicsDBImport and GenotypeGVCFs are run once for each of 18 chromosome
## interval lists. CHROMOSOME_MAP decouples chr1..chr18 labels from the QQ74
## reference's underlying contig names.

joint_genotype() {
    require_file "${SAMPLE_NAME_MAP}"
    require_file "${CHROMOSOME_MAP}"
    [[ "$(chromosome_count)" -eq 18 ]] || {
        echo "ERROR: expected 18 entries in ${CHROMOSOME_MAP}; found $(chromosome_count)" >&2
        exit 1
    }

    local chromosome intervals workspace output_dir
    while IFS=$'\t' read -r chromosome intervals; do
        chromosome="${chromosome%$'\r'}"
        intervals="${intervals%$'\r'}"
        [[ -n "${chromosome}" && "${chromosome}" != \#* ]] || continue
        require_file "${intervals}"

        output_dir="${JOINT_DIR}/chr${chromosome}"
        workspace="${output_dir}/genomicsdb_workspace"
        mkdir -p "${output_dir}"
        [[ ! -e "${workspace}" ]] || {
            echo "ERROR: GenomicsDB workspace already exists: ${workspace}" >&2
            echo "Move or remove it deliberately before rerunning this chromosome." >&2
            exit 1
        }

        gatk GenomicsDBImport \
            --sample-name-map "${SAMPLE_NAME_MAP}" \
            --genomicsdb-workspace-path "${workspace}" \
            --intervals "${intervals}" \
            --reader-threads "${DB_THREADS}"

        gatk GenotypeGVCFs \
            --reference "${REF}" \
            --variant "gendb://${workspace}" \
            --output "${JOINT_DIR}/chr${chromosome}_raw.vcf"
    done < "${CHROMOSOME_MAP}"
}

###############################################################################
## 11. EXTRACT SNPs WITH GATK SELECTVARIANTS
###############################################################################

extract_snps() {
    mkdir -p "${SNP_DIR}"
    local chromosome raw_vcf
    for chromosome in {1..18}; do
        raw_vcf="${JOINT_DIR}/chr${chromosome}_raw.vcf"
        require_file "${raw_vcf}"
        gatk SelectVariants \
            -R "${REF}" \
            -V "${raw_vcf}" \
            --select-type-to-include SNP \
            -O "${SNP_DIR}/chr${chromosome}_snps_only.vcf"
    done
}

###############################################################################
## 12. REMOVE SNPs MISSING MQRankSum OR ReadPosRankSum
###############################################################################

remove_missing_rank_annotations() {
    local chromosome input_vcf output_vcf
    for chromosome in {1..18}; do
        input_vcf="${SNP_DIR}/chr${chromosome}_snps_only.vcf"
        output_vcf="${SNP_DIR}/chr${chromosome}_no_missing_fields.vcf"
        require_file "${input_vcf}"
        awk '$0 ~ /^#/ || ($0 ~ /MQRankSum/ && $0 ~ /ReadPosRankSum/) {print $0}' \
            "${input_vcf}" > "${output_vcf}"
    done
}

###############################################################################
## 13. APPLY GATK HARD FILTERS TO SNPs
###############################################################################
## SNPs are flagged when QD<2.0, FS>60.0, MQ<40.0, MQRankSum<-12.5,
## ReadPosRankSum<-8.0, or SOR>3.0.

hard_filter_snps() {
    local chromosome input_vcf
    for chromosome in {1..18}; do
        input_vcf="${SNP_DIR}/chr${chromosome}_no_missing_fields.vcf"
        require_file "${input_vcf}"
        gatk VariantFiltration \
            -R "${REF}" \
            -V "${input_vcf}" \
            --filter-expression 'QD < 2.0' --filter-name QD_Low \
            --filter-expression 'FS > 60.0' --filter-name FS_High \
            --filter-expression 'MQ < 40.0' --filter-name MQ_Low \
            --filter-expression 'MQRankSum < -12.5' --filter-name MQRankSum_Low \
            --filter-expression 'ReadPosRankSum < -8.0' --filter-name ReadPosRankSum_Low \
            --filter-expression 'SOR > 3.0' --filter-name SOR_High \
            -O "${SNP_DIR}/filtered_chr${chromosome}.vcf"
    done
}

###############################################################################
## 14. RETAIN ONLY PASSING SNPs
###############################################################################
## GATK SelectVariants avoids the INFO-column loss observed in the exploratory
## vcftools recoding step in the original notes.

select_passing_snps() {
    mkdir -p "${FINAL_CHR_DIR}"
    local chromosome input_vcf
    for chromosome in {1..18}; do
        input_vcf="${SNP_DIR}/filtered_chr${chromosome}.vcf"
        require_file "${input_vcf}"
        gatk SelectVariants \
            -R "${REF}" \
            -V "${input_vcf}" \
            --exclude-filtered true \
            -O "${FINAL_CHR_DIR}/pass_chr${chromosome}.vcf"
    done
}

###############################################################################
## 15. REMOVE SNPs WITH MORE THAN 70% MISSING GENOTYPES (BCFTOOLS 1.21)
###############################################################################
## F_MISSING <= 0.7 retains sites with no more than 70% missing data.

filter_missingness() {
    local chromosome input_vcf
    for chromosome in {1..18}; do
        input_vcf="${FINAL_CHR_DIR}/pass_chr${chromosome}.vcf"
        require_file "${input_vcf}"
        bcftools view \
            -i 'F_MISSING<=0.7' \
            -O v \
            -o "${FINAL_CHR_DIR}/chr${chromosome}.vcf" \
            "${input_vcf}"
    done
}

###############################################################################
## 16. MERGE THE 18 CHROMOSOME VCF FILES
###############################################################################

merge_chromosomes() {
    local chromosome inputs=()
    for chromosome in {1..18}; do
        require_file "${FINAL_CHR_DIR}/chr${chromosome}.vcf"
        inputs+=(--INPUT "${FINAL_CHR_DIR}/chr${chromosome}.vcf")
    done

    gatk MergeVcfs \
        "${inputs[@]}" \
        --OUTPUT "${PROJECT_DIR}/combined.filtered.snps.vcf"
}

###############################################################################
## STAGE DISPATCH
###############################################################################

stage="${1:-}"
case "${stage}" in
    download)             download_reads ;;
    trim)                 trim_reads ;;
    reference)            prepare_reference ;;
    align)                align_reads ;;
    read-groups)          add_read_groups_and_mark_duplicates ;;
    index-bam)            index_bams ;;
    call-gvcf)            call_gvcfs ;;
    compress-gvcf)        compress_gvcfs ;;
    sample-map)           create_sample_map ;;
    joint-genotype)       joint_genotype ;;
    extract-snps)         extract_snps ;;
    remove-missing-ranks) remove_missing_rank_annotations ;;
    hard-filter)          hard_filter_snps ;;
    select-pass)          select_passing_snps ;;
    missingness)          filter_missingness ;;
    merge)                merge_chromosomes ;;
    *)
        echo "Usage: $0 {download|trim|reference|align|read-groups|index-bam|call-gvcf|compress-gvcf|sample-map|joint-genotype|extract-snps|remove-missing-ranks|hard-filter|select-pass|missingness|merge}" >&2
        exit 2
        ;;
esac

##Final output file: 91combined.vcf
###############################################################################
## END OF THE ANALYSIS
###############################################################################
