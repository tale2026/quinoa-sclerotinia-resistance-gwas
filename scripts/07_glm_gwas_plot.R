# ==============================================================================
# Script: 05_glm_gwas_plots.R
# ==============================================================================
# Purpose:
#   Import genome-wide association results from the general linear model (GLM)
#   and generate Manhattan and Q-Q plots for resistance to
#   Sclerotinia sclerotiorum in quinoa.
#
# Input:
#   data/gwas/glm_gwas_summary.tsv.gz
#
# Outputs:
#   Manhattan and Q-Q plots for the GLM genome-wide association analysis
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
# 1. Load packages
# ==============================================================================

library(data.table)
library(ggplot2)
library(scales)


# ==============================================================================
# 2. Define input file
# ==============================================================================

glm_gwas_file <- file.path(
  "data",
  "gwas",
  "glm_gwas_summary.tsv.gz"
)


# ==============================================================================
# 3. Confirm that the input file exists
# ==============================================================================

if (!file.exists(glm_gwas_file)) {
  stop(
    "The GLM GWAS file was not found at: ",
    glm_gwas_file,
    "\nRun the script from the repository root and check the filename."
  )
}


# ==============================================================================
# 4. Import GLM GWAS results
# ==============================================================================

glm_gwas <- fread(
  glm_gwas_file,
  sep = "\t",
  header = TRUE,
  na.strings = c("NA", "NaN", ".")
)


# ==============================================================================
# 5. Inspect imported data
# ==============================================================================

dim(glm_gwas)
names(glm_gwas)
str(glm_gwas)
head(glm_gwas, 10)
summary(glm_gwas)

# ==============================================================================
# 6. Manhattan plot
# ==============================================================================

# Keep valid GWAS results
plot_data <- glm_gwas[
  !is.na(chromosome) &
    !is.na(bp) &
    !is.na(p_value) &
    p_value >= 0 &
    p_value <= 1,
  .(
    SNP = locus,
    CHR_NUM = chromosome,
    CHR = as.character(chromosome_label),
    BP = bp,
    P = p_value
  )
]

# Preserve biological chromosome order
chromosome_order <- c(
  paste0("Cq", 1:9, "A"),
  paste0("Cq", 1:9, "B")
)

plot_data[, CHR := factor(
  CHR,
  levels = chromosome_order,
  ordered = TRUE
)]

setorder(plot_data, CHR_NUM, BP)

# Replace exact zero p-values so -log10(P) remains finite
minimum_positive_p <- min(plot_data[P > 0, P], na.rm = TRUE)

plot_data[P == 0, P := minimum_positive_p / 10]

# Calculate chromosome lengths and cumulative offsets
chromosome_info <- plot_data[
  ,
  .(chromosome_length = max(BP, na.rm = TRUE)),
  by = .(CHR_NUM, CHR)
][order(CHR_NUM)]

chromosome_info[, offset :=
                  shift(cumsum(chromosome_length), fill = 0)
]

chromosome_info[, midpoint :=
                  offset + chromosome_length / 2
]

# Add cumulative genomic positions
plot_data <- chromosome_info[
  plot_data,
  on = .(CHR_NUM, CHR)
]

plot_data[, cumulative_BP := BP + offset]
plot_data[, minus_log10_p := -log10(P)]

# Bonferroni threshold based on the number of tested variants
number_of_tests <- nrow(plot_data)
bonferroni_p <- 0.05 / number_of_tests
bonferroni_y <- -log10(bonferroni_p)

cat("Number of tests:", number_of_tests, "\n")
cat("Bonferroni p-value:", format(bonferroni_p, scientific = TRUE), "\n")
cat("Manhattan threshold:", bonferroni_y, "\n")

# Color-blind-friendly Okabe–Ito colors
chromosome_colors <- rep(
  c("#0072B2", "#E69F00"),
  length.out = nrow(chromosome_info)
)

names(chromosome_colors) <- as.character(chromosome_info$CHR)

# Create Manhattan plot
manhattan_plot <- ggplot(
  plot_data,
  aes(
    x = cumulative_BP,
    y = minus_log10_p,
    color = CHR
  )
) +
  geom_point(
    size = 2,
    alpha = 0.75,
    stroke = 0
  ) +
  geom_hline(
    yintercept = bonferroni_y,
    color = "#D55E00",
    linewidth = 0.65,
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = Inf,
    y = bonferroni_y,
    label = paste0(
      "Bonferroni: P = ",
      scientific(bonferroni_p, digits = 2)
    ),
    hjust = 1.03,
    vjust = -0.5,
    size = 3.3,
    color = "black"
  ) +
  scale_color_manual(values = chromosome_colors) +
  scale_x_continuous(
    breaks = chromosome_info$midpoint,
    labels = chromosome_info$CHR,
    expand = expansion(mult = c(0.005, 0.01))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](italic(P))),
    title = ""
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )

manhattan_plot
