# ==============================================================================
# Script: 04_snp_genotype_class_analysis.R
# ==============================================================================
# Purpose:
#   Analyse differences in disease response among genotype classes of
#   significant SNPs identified in the genome-wide association analysis.
#
# Input:
#   data/haplotype_data.xlsx
#   Worksheet: 1
#
# Outputs:
#   Genotype-class frequencies, statistical comparisons, and figures showing
#   disease-response distributions among SNP genotype classes
#
# Author:
#   Swapnil Tale
#
# Last updated:
#   2026-08-28
#
# Note:
#   Run this script from the root directory of the GitHub repository.
# ==============================================================================


# ==============================================================================
# Load packages
# ==============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)


# ==============================================================================
# 1. Define input file and import data
# ==============================================================================

genotype_file <- file.path(
  "data",
  "haplotype_data.xlsx"
)


haplotype_data <- read_excel(
  genotype_file,
  sheet = 1
)


# ------------------------------------------------------------
# 2. Clean genotype columns
# ------------------------------------------------------------

haplotype_data <- haplotype_data %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::starts_with("Cq"),
      ~ {
        cleaned_value <- trimws(as.character(.x))
        
        cleaned_value[
          cleaned_value == "" |
            tolower(cleaned_value) == "missing"
        ] <- NA_character_
        
        cleaned_value
      }
    ),
    audpc = as.numeric(audpc)
  )


# ------------------------------------------------------------
# 3. Convert SNP columns to long format
# ------------------------------------------------------------

snp_long <- haplotype_data %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("Cq"),
    names_to = "SNP",
    values_to = "genotype"
  ) %>%
  dplyr::filter(
    !is.na(genotype),
    !is.na(audpc)
  ) %>%
  tidyr::separate(
    genotype,
    into = c("allele1", "allele2"),
    sep = "/",
    remove = FALSE,
    extra = "merge",
    fill = "right"
  )


# ------------------------------------------------------------
# 4. Count accessions belonging to each genotype
# ------------------------------------------------------------

genotype_counts <- snp_long %>%
  dplyr::count(
    SNP,
    genotype,
    allele1,
    allele2,
    name = "n"
  )


# ------------------------------------------------------------
# 5. Identify the most frequent homozygous genotype
#
# This genotype is treated as the reference genotype.
# ------------------------------------------------------------

reference_genotypes <- genotype_counts %>%
  dplyr::filter(
    allele1 == allele2
  ) %>%
  dplyr::group_by(SNP) %>%
  dplyr::slice_max(
    order_by = n,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    SNP,
    reference_genotype = genotype
  )


# ------------------------------------------------------------
# 6. Assign genotype classes
# ------------------------------------------------------------

genotype_counts <- genotype_counts %>%
  dplyr::left_join(
    reference_genotypes,
    by = "SNP"
  ) %>%
  dplyr::mutate(
    genotype_class = dplyr::case_when(
      allele1 != allele2 ~
        "Heterozygous",
      
      genotype == reference_genotype ~
        "Homozygous reference",
      
      TRUE ~
        "Homozygous alternate"
    ),
    
    genotype_class = factor(
      genotype_class,
      levels = c(
        "Homozygous reference",
        "Heterozygous",
        "Homozygous alternate"
      )
    ),
    
    group_label = paste0(
      genotype,
      "\nn = ",
      n
    )
  )


# ------------------------------------------------------------
# 7. Add genotype classes to the plotting data
# ------------------------------------------------------------

plot_data <- snp_long %>%
  dplyr::left_join(
    genotype_counts %>%
      dplyr::select(
        SNP,
        genotype,
        genotype_class,
        group_label
      ),
    by = c("SNP", "genotype")
  )


# ------------------------------------------------------------
# 8. Define SNP order
# ------------------------------------------------------------

snp_order <- c(
  "Cq5B_58690487",
  "Cq5B_58849261",
  "Cq5B_59088400",
  "Cq5B_59346122",
  "Cq5B_59716352",
  "Cq9A_51113865"
)

plot_data <- plot_data %>%
  dplyr::mutate(
    SNP = factor(
      SNP,
      levels = snp_order
    )
  )

genotype_counts <- genotype_counts %>%
  dplyr::mutate(
    SNP = factor(
      SNP,
      levels = snp_order
    )
  )


# ------------------------------------------------------------
# 9. Define panel labels
# ------------------------------------------------------------

panel_labels <- c(
  "Cq5B_58690487" =
    "bold(c)~plain('Cq5B_58690487')",
  
  "Cq5B_58849261" =
    "bold(d)~plain('Cq5B_58849261')",
  
  "Cq5B_59088400" =
    "bold(e)~plain('Cq5B_59088400')",
  
  "Cq5B_59346122" =
    "bold(f)~plain('Cq5B_59346122')",
  
  "Cq5B_59716352" =
    "bold(g)~plain('Cq5B_59716352')",
  
  "Cq9A_51113865" =
    "bold(h)~plain('Cq9A_51113865')"
)


# ------------------------------------------------------------
# 10. Position genotype and sample-size labels
# ------------------------------------------------------------

label_y <- max(
  plot_data$audpc,
  na.rm = TRUE
) + 2

genotype_label_data <- genotype_counts %>%
  dplyr::mutate(
    label_y = label_y
  )


# ------------------------------------------------------------
# 11. Set upper y-axis limit
# ------------------------------------------------------------

plot_upper_limit <- label_y + 6


# ------------------------------------------------------------
# 12. Create plot
# ------------------------------------------------------------

haplotype_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = genotype_class,
    y = audpc,
    fill = genotype_class
  )
) +
  
  # Boxplots
  ggplot2::geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    alpha = 0.80,
    colour = "black",
    linewidth = 0.35
  ) +
  
  # Individual biological replicates
  ggplot2::geom_jitter(
    width = 0.10,
    height = 0,
    shape = 21,
    size = 1.10,
    stroke = 0.20,
    colour = "black",
    fill = "grey35",
    alpha = 0.70
  ) +
  
  # Genotype and sample-size labels
  ggplot2::geom_text(
    data = genotype_label_data,
    ggplot2::aes(
      x = genotype_class,
      y = label_y,
      label = group_label
    ),
    inherit.aes = FALSE,
    family = "Arial",
    size = 2.5,
    lineheight = 0.95,
    vjust = 0
  ) +
  
  # Six panels arranged in two columns
  ggplot2::facet_wrap(
    ~ SNP,
    ncol = 2,
    drop = TRUE,
    labeller = ggplot2::labeller(
      SNP = ggplot2::as_labeller(
        panel_labels,
        default = ggplot2::label_parsed
      )
    )
  ) +
  
  # Colour-blind-friendly palette
  ggplot2::scale_fill_manual(
    values = c(
      "Homozygous reference" = "#0072B2",
      "Heterozygous"         = "#E69F00",
      "Homozygous alternate" = "#CC79A7"
    ),
    breaks = c(
      "Homozygous reference",
      "Heterozygous",
      "Homozygous alternate"
    ),
    drop = FALSE,
    name = "Genotype class"
  ) +
  
  ggplot2::scale_x_discrete(
    drop = FALSE
  ) +
  
  ggplot2::scale_y_continuous(
    limits = c(0, plot_upper_limit),
    breaks = seq(
      0,
      45,
      by = 5
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.01)
    )
  ) +
  
  ggplot2::labs(
    x = "Genotype class",
    y = "AUDPC"
  ) +
  
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      title.position = "left",
      nrow = 1,
      byrow = TRUE
    )
  ) +
  
  ggplot2::theme_bw(
    base_size = 9,
    base_family = "Arial"
  ) +
  
  ggplot2::theme(
    # Facet labels
    strip.background = ggplot2::element_blank(),
    
    strip.text = ggplot2::element_text(
      face = "plain",
      size = 8.5,
      colour = "black",
      margin = ggplot2::margin(
        t = 2,
        b = 2
      )
    ),
    
    # Axis titles
    axis.title.x = ggplot2::element_text(
      size = 9,
      margin = ggplot2::margin(t = 5)
    ),
    
    axis.title.y = ggplot2::element_text(
      size = 9,
      margin = ggplot2::margin(r = 5)
    ),
    
    # Axis labels
    axis.text.x = ggplot2::element_text(
      angle = 40,
      size = 7,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    
    axis.text.y = ggplot2::element_text(
      size = 7.5,
      colour = "black"
    ),
    
    axis.ticks = ggplot2::element_line(
      linewidth = 0.30,
      colour = "black"
    ),
    
    # Unboxed panels
    panel.border = ggplot2::element_blank(),
    
    axis.line = ggplot2::element_line(
      colour = "black",
      linewidth = 0.45
    ),
    
    panel.grid.major.x = ggplot2::element_blank(),
    
    panel.grid.minor = ggplot2::element_blank(),
    
    panel.grid.major.y = ggplot2::element_line(
      colour = "grey85",
      linewidth = 0.30
    ),
    
    # Legend
    legend.position = "bottom",
    
    legend.title = ggplot2::element_text(
      size = 8,
      face = "bold"
    ),
    
    legend.text = ggplot2::element_text(
      size = 8
    ),
    
    legend.key.size = grid::unit(
      4,
      "mm"
    ),
    
    # Panel spacing
    panel.spacing.x = grid::unit(
      6,
      "mm"
    ),
    
    panel.spacing.y = grid::unit(
      5,
      "mm"
    ),
    
    plot.margin = ggplot2::margin(
      t = 4,
      r = 4,
      b = 4,
      l = 4,
      unit = "mm"
    )
  )


# ------------------------------------------------------------
# 13. Display plot
# ------------------------------------------------------------

haplotype_plot


# ------------------------------------------------------------
# 14. Save files
# ------------------------------------------------------------

# TIFF: 600 dpi
ggplot2::ggsave(
  filename = file.path(
    output_directory,
    "Haplotype_AUDPC_boxplots3.tiff"
  ),
  plot = haplotype_plot,
  device = "tiff",
  width = 130,
  height = 190,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# ------------------------------------------------------------
# END OF THE ANALYSIS
# ------------------------------------------------------------
