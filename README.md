DATA AND ANALYSIS SCRIPTS SUPPORTING A GENOME-WIDE ASSOCIATION STUDY
OF RESISTANCE TO SCLEROTINIA SCLEROTIORUM IN QUINOA
===============================================================================

Version: 1.0.0
Release date: 2026-08-28
Zenodo DOI: https://doi.org/10.5281/zenodo.22146648


1. ASSOCIATED MANUSCRIPT
-------------------------------------------------------------------------------

Title:
Genetic Dissection of Sclerotinia Stem Rot Resistance in Quinoa Reveals Resistance-Associated Loci and Candidate Defense Genes

Authors:
Swapnil Tale, Severin Einspanier, Remco Stam, and Nazgol Emrani


2. DATASET DESCRIPTION
-------------------------------------------------------------------------------

This record contains phenotypic data and analysis scripts supporting
a study of resistance to Sclerotinia sclerotiorum in a diversity panel of
91 quinoa accessions.

The deposited files include data from detached-leaf and whole-plant
stem-inoculation experiments, lesion-length measurements, area under the
disease progress curve (AUDPC) values, and scripts used for statistical analysis and visualization.


3. SOURCE DATA AND PROVENANCE
-------------------------------------------------------------------------------

Raw sequencing reads:
The genotype data were generated from publicly available sequencing reads
associated with NCBI BioProject PRJNA673789:

https://www.ncbi.nlm.nih.gov/bioproject/PRJNA673789

Reference genome:
The QQ74 v2 reference genome used for read alignment and variant calling is
available through the CoGe database under genome ID 60716:

https://genomevolution.org/coge/GenomeInfo.pl?gid=60716

Reference-genome publication:
Rey, E., Maughan, P. J., Maumus, F., Lewis, D., Wilson, L., Fuller, J., Schmöckel, S. M., Jellen, E. N., Tester, M., & Jarvis, D. E. (2023). A chromosome-scale assembly of the quinoa genome provides insights into the structure and dynamics of its subgenomes. Communications Biology, 6(1), 1263.

The raw sequencing reads and reference-genome files are not duplicated in
the Zenodo record. They remain available through their original repositories.


4. FILE DESCRIPTIONS
-------------------------------------------------------------------------------

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


5. ANALYSIS SCRIPTS
-------------------------------------------------------------------------------

01_detached_leaf_assay.R
    Imports and prepares detached-leaf assay data and performs the associated
    descriptive and statistical analyses.

02_stem_lesion_age_test.R
    Analyses the preliminary stem-lesion experiment comparing disease responses
    in six- and nine-week-old quinoa plants.

03_stem_lesion_audpc_analysis.R
    Analyses lesion length at 7, 12, and 17 days post-inoculation and area under
    the disease progress curve (AUDPC) among the 91 quinoa accessions.

04_variant_calling.sh
    Performs read processing, alignment, and variant-calling steps used to
    generate the genotype VCF from the publicly available sequencing data.
    See the script header for software requirements and input paths.

05_glm_gwas_hail.ipynb
    Jupyter notebook implementing genotype filtering and the general linear
    model genome-wide association analysis using Hail.

06_mlm_gwas_gcta.sh
    Prepares the Hail-exported PLINK data for GCTA, constructs the genomic
    relationship matrix, calculates principal components, and performs the
    mixed linear model association analysis using PC1 as a covariate.

07_glm_gwas_plot.R
    Imports the GLM association summary statistics and generates the
    corresponding Manhattan plot.

08_mlm_gwas_plot.R
    Imports the MLM association summary statistics and generates the
    corresponding Manhattan plot

09_snp_genotype_class_analysis.R
    Analyses phenotypic differences among genotype classes of SNPs identified
    in the genome-wide association analyses and generates the associated
    genotype-class figures.

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

1. Download the Zenodo record and the associated data files.
2. Consult README.txt for the required input-file locations.
3. Run 04_variant_calling.sh to reproduce the genotype VCF from the publicly
   available sequencing reads, if complete regeneration is required.
4. Run 05_glm_gwas_hail.ipynb to filter the genotype data and conduct the GLM
   association analysis.
5. Run 06_mlm_gwas_gcta.sh to construct the genomic relationship matrix,
   calculate principal components, and conduct the MLM association analysis.
6. Run 01_detached_leaf_assay.R, 02_stem_lesion_age_test.R, and
   03_stem_lesion_audpc_analysis.R for the phenotypic analyses.
7. Run 07_glm_gwas_plot.R and 08_mlm_gwas_plot.R to reproduce the Manhattan
   and quantile-quantile plots.
8. Run 09_snp_genotype_class_analysis.R to reproduce the SNP genotype-class
   comparisons and figures.

The R scripts should be run from the repository root directory. The Bash
scripts require a Linux environment and the software specified in their
headers. The Hail analysis requires a Python environment with Hail and
Jupyter installed.


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

Tale S, Einspanier S, Stam R, Emrani N (2026). Data and analysis scripts
supporting a genome-wide association study of resistance to Sclerotinia
sclerotiorum in quinoa. Zenodo. https://doi.org/10.5281/zenodo.22146648 


11. CONTACT
-------------------------------------------------------------------------------

Swapnil Tale
Crop Genetics, University of Rostock, Justus-von-Liebig-Weg 8, D-18059 Rostock, Germany 
Email: swapnil.tale@uni-rostock.de
ORCID: https://orcid.org/0009-0007-5960-7084
