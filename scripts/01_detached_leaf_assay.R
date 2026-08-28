# ==============================================================================
# Script: 01_detached_leaf_assay.R
# ==============================================================================
# Purpose:
#   Analyse detached leaf assay data used to evaluate the pathogenicity of
#   Sclerotinia sclerotiorum and differences in lesion area among quinoa
#   accessions.
#
# Input:
#   data/phenotypic_data.xlsx
#   Worksheet: Detached_leaf
#
# Outputs:
#   Statistical results and figures from the detached leaf assay
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

#------------------------------------------------------------
# 1. Load packages
#------------------------------------------------------------

library(readxl)
library(dplyr)
library(ggplot2)
library(broom)
library(car)
library(emmeans)
library(multcomp)
library(multcompView)
library(openxlsx)


#------------------------------------------------------------
# 2. Import and prepare data
#------------------------------------------------------------
library(readxl)

phenotypic_file <- file.path(
  "data",
  "phenotypic_data.xlsx"
)

detached_leaf_data <- read_excel(
  phenotypic_file,
  sheet = 1
) %>%
  dplyr::select(
    accession,
    Lesion_Size
  ) %>%
  dplyr::filter(
    !is.na(accession),
    !is.na(Lesion_Size)
  ) %>%
  dplyr::mutate(
    accession = as.character(accession),
    Lesion_Size = as.numeric(Lesion_Size)
  )

# Preserve the order in the Excel file
accession_order <- unique(detached_leaf_data$accession)

detached_leaf_data <- detached_leaf_data %>%
  dplyr::mutate(
    accession = factor(
      accession,
      levels = accession_order
    )
  )

str(detached_leaf_data)


#------------------------------------------------------------
# 3. Descriptive statistics
#------------------------------------------------------------

descriptive_results <- detached_leaf_data %>%
  dplyr::group_by(accession) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean = mean(Lesion_Size),
    SD = sd(Lesion_Size),
    SE = SD / sqrt(n),
    median = median(Lesion_Size),
    minimum = min(Lesion_Size),
    maximum = max(Lesion_Size),
    .groups = "drop"
  )


#------------------------------------------------------------
# 4. One-way ANOVA
#------------------------------------------------------------

anova_model <- aov(
  Lesion_Size ~ accession,
  data = detached_leaf_data
)

anova_results <- broom::tidy(anova_model)

print(summary(anova_model))


#------------------------------------------------------------
# 5. Assumption tests
#------------------------------------------------------------

# Normality should be tested on model residuals
shapiro_results <- broom::tidy(
  shapiro.test(
    residuals(anova_model)
  )
)

# Homogeneity of variance
levene_results <- broom::tidy(
  car::leveneTest(
    Lesion_Size ~ accession,
    data = detached_leaf_data,
    center = median
  )
)

shapiro_results
levene_results


#------------------------------------------------------------
# 6. Tukey-adjusted pairwise comparisons
#------------------------------------------------------------

emm <- emmeans::emmeans(
  anova_model,
  ~ accession
)

tukey_results <- emmeans::contrast(
  emm,
  method = "pairwise",
  adjust = "tukey"
) %>%
  as.data.frame()

tukey_results

# Display only significant comparisons
significant_tukey_results <- tukey_results %>%
  dplyr::filter(p.value < 0.05) %>%
  dplyr::arrange(p.value)

significant_tukey_results


#------------------------------------------------------------
# 7. Tukey compact-letter display
#------------------------------------------------------------

letters_df <- multcomp::cld(
  emm,
  alpha = 0.05,
  adjust = "tukey",
  Letters = base::letters,
  sort = FALSE
) %>%
  as.data.frame() %>%
  dplyr::mutate(
    .group = trimws(.group)
  )

# The use of base::letters is important because an existing
# object named "letters" is present in the R environment.


#------------------------------------------------------------
# 8. Calculate positions for significance letters
#------------------------------------------------------------

label_offset <- 0.08 * diff(
  range(
    detached_leaf_data$Lesion_Size,
    na.rm = TRUE
  )
)

label_positions <- detached_leaf_data %>%
  dplyr::group_by(accession) %>%
  dplyr::summarise(
    y_pos = max(
      Lesion_Size,
      na.rm = TRUE
    ) + label_offset,
    .groups = "drop"
  )

letters_df <- letters_df %>%
  dplyr::left_join(
    label_positions,
    by = "accession"
  ) %>%
  dplyr::arrange(accession)

# Inspect the final means and letters
letters_df %>%
  dplyr::select(
    accession,
    emmean,
    SE,
    .group,
    y_pos
  )


#------------------------------------------------------------
# 9. Save statistical results to Excel
#------------------------------------------------------------

statistics_file <- file.path(
  output_directory,
  "Detached_Leaf_Assay_Statistics.xlsx"
)

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Descriptive")
openxlsx::writeData(
  wb,
  "Descriptive",
  descriptive_results
)

openxlsx::addWorksheet(wb, "Shapiro_residuals")
openxlsx::writeData(
  wb,
  "Shapiro_residuals",
  shapiro_results
)

openxlsx::addWorksheet(wb, "Levene")
openxlsx::writeData(
  wb,
  "Levene",
  levene_results
)

openxlsx::addWorksheet(wb, "ANOVA")
openxlsx::writeData(
  wb,
  "ANOVA",
  anova_results
)

openxlsx::addWorksheet(wb, "Tukey_HSD")
openxlsx::writeData(
  wb,
  "Tukey_HSD",
  tukey_results
)

openxlsx::addWorksheet(wb, "Tukey_significant")
openxlsx::writeData(
  wb,
  "Tukey_significant",
  significant_tukey_results
)

openxlsx::addWorksheet(wb, "Tukey_letters")
openxlsx::writeData(
  wb,
  "Tukey_letters",
  letters_df
)

openxlsx::saveWorkbook(
  wb,
  file = statistics_file,
  overwrite = TRUE
)


#------------------------------------------------------------
# 10. Journal-ready plot
#------------------------------------------------------------

p_detached_leaf <- ggplot(
  detached_leaf_data,
  aes(
    x = accession,
    y = Lesion_Size
  )
) +
  
  # Individual biological observations
  geom_point(
    position = position_jitter(
      width = 0.10,
      height = 0,
      seed = 123
    ),
    shape = 21,
    size = 1.8,
    stroke = 0.35,
    colour = "black",
    fill = "grey70",
    alpha = 0.85
  ) +
  
  # Mean ± standard error
  stat_summary(
    fun.data = ggplot2::mean_se,
    geom = "errorbar",
    width = 0.16,
    colour = "black",
    linewidth = 0.5
  ) +
  
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 3,
    colour = "black"
  ) +
  
  # Tukey significance letters
  geom_text(
    data = letters_df,
    aes(
      x = accession,
      y = y_pos,
      label = .group
    ),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "bold",
    size = 3
  ) +
  
  labs(
    x = "Accession",
    y = expression("Lesion area (cm"^2*")")
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0.02, 0.14)
    )
  ) +
  
  theme_classic(
    base_size = 9,
    base_family = "Arial"
  ) +
  
  theme(
    axis.title = element_text(
      colour = "black",
      size = 9
    ),
    axis.text = element_text(
      colour = "black",
      size = 8
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.4
    ),
    plot.margin = margin(
      4, 4, 4, 4,
      unit = "mm"
    )
  )

p_detached_leaf


#------------------------------------------------------------
# 11. Save figure as 600-dpi TIFF
#------------------------------------------------------------

ggsave(
  filename = file.path(
    output_directory,
    "Detached_leaf_lesion_area.tiff"
  ),
  plot = p_detached_leaf,
  device = "tiff",
  width = 84,
  height = 82,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)
#------------------------------------------------------------
#  END OF THE ANALYSIS
#------------------------------------------------------------