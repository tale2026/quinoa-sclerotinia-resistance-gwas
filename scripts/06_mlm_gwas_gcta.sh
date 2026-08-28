#!/usr/bin/env bash

# ==============================================================================
# Script: 06_mlm_gwas_gcta.sh
# ==============================================================================
# Purpose:
#   Perform a mixed linear model genome-wide association analysis of mean AUDPC
#   values using GCTA and the high-confidence quinoa SNP dataset exported from
#   Hail in PLINK binary format.
#
# Genotype input:
#   data/genotype/plink_91_accessions_filtered/
#     quinoa_91_accessions_filtered.bed
#     quinoa_91_accessions_filtered.bim
#     quinoa_91_accessions_filtered.fam
#
# Phenotype input:
#   data/phenotype/mean_audpc.txt
#
# Analysis:
#   1. Convert QQ74 contig identifiers to chromosome numbers 1-18.
#   2. Construct a genomic relationship matrix.
#   3. Calculate the first 10 principal components.
#   4. Use PC1 as a quantitative covariate.
#   5. Perform the MLM association analysis using GCTA-MLMA.
#
# Outputs:
#   Modified PLINK files, genomic relationship matrix, principal components,
#   PC1 covariate file, and MLM GWAS results
#
# Author:
#   Swapnil Tale
#
# Last updated:
#   2026-08-28
#
# Requirements:
#   Bash and GCTA
# ==============================================================================

set -euo pipefail


# ==============================================================================
# 1. Determine repository location
# ==============================================================================

SCRIPT_DIR="$(
    cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &&
    pwd
)"

REPO_ROOT="$(
    cd -- "${SCRIPT_DIR}/.." &&
    pwd
)"


# ==============================================================================
# 2. Configure input and output paths
# ==============================================================================

GCTA="${GCTA:-gcta64}"

GENOTYPE_DIR="${GENOTYPE_DIR:-${REPO_ROOT}/data/genotype/plink_91_accessions_filtered}"

INPUT_BASENAME="${INPUT_BASENAME:-quinoa_91_accessions_filtered}"

INPUT_PREFIX="${INPUT_PREFIX:-${GENOTYPE_DIR}/${INPUT_BASENAME}}"

PHENOTYPE_INPUT="${PHENOTYPE_INPUT:-${REPO_ROOT}/data/phenotype/mean_audpc_gcta.tsv}"

OUTPUT_DIR="${OUTPUT_DIR:-${REPO_ROOT}/results/mlm_gwas_gcta}"

GCTA_PREFIX="${GCTA_PREFIX:-${OUTPUT_DIR}/quinoa_91_accessions_gcta}"

PHENOTYPE_GCTA="${PHENOTYPE_GCTA:-${OUTPUT_DIR}/mean_audpc_noheader.txt}"

GRM_PREFIX="${GRM_PREFIX:-${OUTPUT_DIR}/quinoa_91_accessions_grm}"

PCA_PREFIX="${PCA_PREFIX:-${OUTPUT_DIR}/quinoa_91_accessions_pcs}"

PC1_FILE="${PC1_FILE:-${OUTPUT_DIR}/pc1.qcovar}"

RESULT_PREFIX="${RESULT_PREFIX:-${OUTPUT_DIR}/mlm_gwas_pc1}"

mkdir -p "${OUTPUT_DIR}"


# ==============================================================================
# 3. Check required programs and files
# ==============================================================================

require_file() {
    [[ -s "$1" ]] || {
        echo "ERROR: required file is absent or empty: $1" >&2
        exit 1
    }
}

command -v "${GCTA}" >/dev/null 2>&1 || {
    echo "ERROR: GCTA executable not found: ${GCTA}" >&2
    exit 1
}

require_file "${INPUT_PREFIX}.bed"
require_file "${INPUT_PREFIX}.bim"
require_file "${INPUT_PREFIX}.fam"
require_file "${PHENOTYPE_INPUT}"


# ==============================================================================
# 4. Prepare the Hail-exported PLINK dataset for GCTA
# ==============================================================================
# The BED genotype data and FAM sample information are retained unchanged.
# QQ74 contig identifiers in the BIM file are subsequently converted to
# chromosome numbers 1-18.

cp "${INPUT_PREFIX}.bed" "${GCTA_PREFIX}.bed"
cp "${INPUT_PREFIX}.fam" "${GCTA_PREFIX}.fam"

awk 'BEGIN {OFS="\t"}
{
    if      ($1 == "lcl|Contig2_pilon")   $1 = 1;
    else if ($1 == "lcl|Contig25_pilon")  $1 = 2;
    else if ($1 == "lcl|Contig37_pilon")  $1 = 3;
    else if ($1 == "lcl|Contig4_pilon")   $1 = 4;
    else if ($1 == "lcl|Contig18_pilon")  $1 = 5;
    else if ($1 == "lcl|Contig19_pilon")  $1 = 6;
    else if ($1 == "lcl|Contig20_pilon")  $1 = 7;
    else if ($1 == "lcl|Contig133_pilon") $1 = 8;
    else if ($1 == "lcl|Contig41_pilon")  $1 = 9;
    else if ($1 == "lcl|Contig1_pilon")   $1 = 10;
    else if ($1 == "lcl|Contig27_pilon")  $1 = 11;
    else if ($1 == "lcl|Contig43_pilon")  $1 = 12;
    else if ($1 == "lcl|Contig32_pilon")  $1 = 13;
    else if ($1 == "lcl|Contig12_pilon")  $1 = 14;
    else if ($1 == "lcl|Contig92_pilon")  $1 = 15;
    else if ($1 == "lcl|Contig28_pilon")  $1 = 16;
    else if ($1 == "lcl|Contig51_pilon")  $1 = 17;
    else if ($1 == "lcl|Contig48_pilon")  $1 = 18;
    else {
        print "ERROR: unmapped reference contig: " $1 > "/dev/stderr";
        exit 1;
    }

    print;
}' "${INPUT_PREFIX}.bim" > "${GCTA_PREFIX}.bim"

## Confirm that all output chromosome identifiers are numerical values 1-18.
awk '
    $1 !~ /^[0-9]+$/ || $1 < 1 || $1 > 18 {
        print "ERROR: invalid chromosome identifier in modified BIM: " $1 \
            > "/dev/stderr";
        exit 1;
    }
' "${GCTA_PREFIX}.bim"

###############################################################################
## 5. FORMAT THE PHENOTYPE FILE FOR GCTA
###############################################################################
## The original workflow assigned FID=0 and removed the header. This is valid
## only when column 1 of the HAIL-exported FAM file also contains FID=0.
## Output format: FID, IID, mean_AUDPC; no header.

awk 'BEGIN {OFS="\t"}
     NR == 1 {next}
     NF >= 3 {$1 = 0; print $1, $2, $3}' \
    "${PHENOTYPE_INPUT}" > "${PHENOTYPE_GCTA}"

require_file "${PHENOTYPE_GCTA}"

###############################################################################
## 6. VERIFY PHENOTYPE AND GENOTYPE SAMPLE IDENTIFIERS
###############################################################################
## Every FID/IID combination in the phenotype file must occur in the FAM file.

awk 'BEGIN {OFS="\t"} {print $1, $2}' "${GCTA_PREFIX}.fam" |
    sort -u > "${RESULT_PREFIX}.fam_ids.tmp"

awk 'BEGIN {OFS="\t"} {print $1, $2}' "${PHENOTYPE_GCTA}" |
    sort -u > "${RESULT_PREFIX}.phenotype_ids.tmp"

if ! comm -23 \
    "${RESULT_PREFIX}.phenotype_ids.tmp" \
    "${RESULT_PREFIX}.fam_ids.tmp" \
    > "${RESULT_PREFIX}.unmatched_ids.tmp"; then
    echo "ERROR: sample-identifier comparison failed." >&2
    exit 1
fi

if [[ -s "${RESULT_PREFIX}.unmatched_ids.tmp" ]]; then
    echo "ERROR: phenotype FID/IID values not found in ${GCTA_PREFIX}.fam:" >&2
    cat "${RESULT_PREFIX}.unmatched_ids.tmp" >&2
    echo "Check whether the FAM file uses FID=0 before running this analysis." >&2
    exit 1
fi

rm -f \
    "${RESULT_PREFIX}.fam_ids.tmp" \
    "${RESULT_PREFIX}.phenotype_ids.tmp" \
    "${RESULT_PREFIX}.unmatched_ids.tmp"

###############################################################################
## 7. CONSTRUCT THE GENOMIC RELATIONSHIP MATRIX
###############################################################################

"${GCTA}" \
    --bfile "${GCTA_PREFIX}" \
    --make-grm \
    --out "${GRM_PREFIX}"

###############################################################################
##8. CALCULATE THE FIRST 10 PRINCIPAL COMPONENTS
###############################################################################

"${GCTA}" \
    --grm "${GRM_PREFIX}" \
    --pca 10 \
    --out "${PCA_PREFIX}"

###############################################################################
## 9. PREPARE PC1 AS A QUANTITATIVE COVARIATE
###############################################################################
## GCTA eigenvector output contains FID, IID, PC1, PC2, ..., PC10.
## Only PC1 was included in the final MLM reported by the original workflow.

require_file "${PCA_PREFIX}.eigenvec"

awk 'BEGIN {OFS="\t"} {print $1, $2, $3}' \
    "${PCA_PREFIX}.eigenvec" > "${PC1_FILE}"

###############################################################################
## 10. PERFORM THE MIXED LINEAR MODEL GWAS
###############################################################################
## The model includes the genomic relationship matrix as the random-effect
## covariance structure and PC1 as a quantitative covariate.

"${GCTA}" \
    --bfile "${GCTA_PREFIX}" \
    --grm "${GRM_PREFIX}" \
    --pheno "${PHENOTYPE_GCTA}" \
    --qcovar "${PC1_FILE}" \
    --mlma \
    --out "${RESULT_PREFIX}"

echo "GCTA MLM GWAS completed: ${RESULT_PREFIX}.mlma"
###############################################################################
## END OF THE ANALYSIS
###############################################################################
