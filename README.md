DATA AND ANALYSIS SCRIPTS SUPPORTING A GENOME-WIDE ASSOCIATION STUDY
OF RESISTANCE TO SCLEROTINIA SCLEROTIORUM IN QUINOA
===============================================================================

Version: 1.0.0
Release date: 2026-08-28
Zenodo DOI: https://doi.org/10.5281/zenodo.22146648


1. ASSOCIATED MANUSCRIPT
-------------------------------------------------------------------------------

Title:
Genome-Wide Association Study of Resistance to Sclerotinia sclerotiorum
in Quinoa Identifies Resistance-Associated Loci and Defense-Related
Candidate Genes

Authors:
Swapnil Tale, Severin Einspanier, Remco Stam, and Nazgol Emrani


2. DATASET DESCRIPTION
-------------------------------------------------------------------------------

This record contains phenotypic data, derived genotype data, genome-wide
association study (GWAS) summary statistics, and analysis scripts supporting
a study of resistance to Sclerotinia sclerotiorum in a diversity panel of
91 quinoa accessions.

The deposited files include data from detached-leaf and whole-plant
stem-inoculation experiments, lesion-length measurements, area under the
disease progress curve (AUDPC) values, filtered variant data, PLINK binary
files used for mixed linear model analysis, GLM and MLM association summary
statistics, and scripts used for statistical analysis and visualization.


3. SOURCE DATA AND PROVENANCE
-------------------------------------------------------------------------------

Raw sequencing reads:
The genotype data were generated from publicly available sequencing reads
associated with NCBI BioProject PRJNA673789:

https://www.ncbi.nlm.nih.gov/bioproject/PRJNA673789

Source publication:
Patiranage et al. (2022)

Reference genome:
The QQ74 v2 reference genome used for read alignment and variant calling is
available through the CoGe database under genome ID 60716:

https://genomevolution.org/coge/GenomeInfo.pl?gid=60716

Reference-genome publication:
Rey et al. (2023)

The raw sequencing reads and reference-genome files are not duplicated in
this Zenodo record. They remain available through their original repositories.


4. FILE DESCRIPTIONS
-------------------------------------------------------------------------------

README.txt
    Description of the dataset, file contents, provenance, and analysis
    workflow.

phenotypic_data.xlsx
    Processed phenotypic data used in the statistical analyses.

    Worksheet 1:
    Detached-leaf assay data.

    Worksheet 2:
    Preliminary stem-lesion experiment comparing disease responses in
    six- and nine-week-old quinoa plants.

    Worksheet 3:
    Stem-lesion measurements at 7, 12, and 17 days post-inoculation,
    AUDPC values, and accession-origin information for the 91-accession
    diversity panel.

haplotype_data.xlsx
    Phenotypic values and genotype classes for SNP genotype-class analyses.
    Missing genotype calls are represented as NA.

quinoa_91_accessions.vcf[.gz or .bgz]
    Derived variant-call data for the 91 quinoa accessions.

quinoa_91_accessions_filtered.bed
quinoa_91_accessions_filtered.bim
quinoa_91_accessions_filtered.fam
    Filtered PLINK binary genotype dataset used to construct the genomic
    relationship matrix, calculate principal components, and perform the
    GCTA mixed linear model association analysis.

glm_gwas_summary.tsv.gz
    Compressed tab-separated GLM genome-wide association summary statistics.

mlm_gwas_summary.tsv.gz
    Compressed tab-separated MLM genome-wide association summary statistics.

analysis_scripts_v1.0.0.zip
    Fixed version of the R, Python, and Bash analysis scripts associated with
    this dataset. The actively maintained scripts are available through the
    GitHub repository listed above.

Supplementary_Tables.xlsx
    Table S1 provides the country of origin and geographic region of the
    91 quinoa accessions. Table S2 lists genes located within 35 kb upstream
    and downstream of the SNPs identified using the mixed linear model.

Supplementary_Data_S1.xlsx
    Detailed statistical results from the preliminary experiment comparing
    disease responses in six- and nine-week-old plants.

Supplementary_Data_S2.xlsx
    ANOVA, Kruskal-Wallis, Tukey HSD, and Dunn post hoc test results for
    lesion length at 7, 12, and 17 days post-inoculation and for AUDPC.


5. ANALYSIS SCRIPTS
-------------------------------------------------------------------------------

01_detached_leaf_assay.R
    Statistical analysis of the detached-leaf assay.

02_stem_lesion_age_test.R
    Analysis of the preliminary experiment comparing six- and
    nine-week-old plants.

03_stem_lesion_audpc_analysis.R
    Analysis of lesion-length progression and AUDPC among the 91 accessions.

04_snp_genotype_class_analysis.R
    Analysis of disease-response differences among SNP genotype classes.

05_glm_gwas_plots.R
    Import, processing, and visualization of GLM GWAS results.

06_mlm_gwas_gcta.sh
    Construction of the genomic relationship matrix, principal-component
    analysis, and MLM association analysis using GCTA.

[INSERT HAIL SCRIPT NAME].py
    Filtering of genotype data and export of the PLINK BED, BIM, and FAM
    files used by GCTA.


6. SOFTWARE
-------------------------------------------------------------------------------

Phenotypic analyses and visualizations were conducted in R version 4.4.2.
The analyses used packages including lme4, boot, car, emmeans, agricolae,
multcompView, readxl, dplyr, tidyr, data.table, broom, openxlsx, ggplot2,
scales, and GGally.

Genotype processing was conducted using Hail. Mixed linear model association
analysis was conducted using GCTA.

Exact package and software versions are documented in the analysis scripts
and session-information file, where provided.


7. REPRODUCING THE ANALYSES
-------------------------------------------------------------------------------

1. Download this Zenodo record.
2. Download or clone the associated GitHub repository.
3. Place the downloaded data files in the repository directories described
   in the GitHub README.
4. Run the R scripts from the root directory of the repository.
5. Run the Hail genotype-processing script to reproduce the filtered PLINK
   files, if required.
6. Run 06_mlm_gwas_gcta.sh in a Linux environment with GCTA available in PATH.
7. Consult the individual script headers for input files, outputs, and
   additional software requirements.

The deposited PLINK files allow the GCTA analysis to be reproduced without
rerunning the complete upstream genotype-processing workflow.


8. MISSING DATA
-------------------------------------------------------------------------------

Missing phenotypic or genotype observations are represented as NA or empty
cells, as documented in the relevant files. Analysis scripts explicitly
exclude missing observations where required.


9. LICENCE
-------------------------------------------------------------------------------

The data are distributed under the Creative Commons Attribution 4.0
International licence (CC BY 4.0), unless otherwise indicated.

The publicly available source sequencing reads and reference-genome files
remain subject to the terms of their original repositories.


10. CITATION
-------------------------------------------------------------------------------

Please cite the specific Zenodo version used:

Tale S, Einspanier S, Stam R, Emrani N (2026). Data and analysis scripts
supporting a genome-wide association study of resistance to Sclerotinia
sclerotiorum in quinoa. Zenodo. https://doi.org/[INSERT DOI]


11. CONTACT
-------------------------------------------------------------------------------

Swapnil Tale
[INSERT AFFILIATION]
Email: [INSERT EMAIL ADDRESS]
ORCID: [INSERT ORCID]
