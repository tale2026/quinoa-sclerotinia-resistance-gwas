# ==============================================================================
# Script: 02_stem_lesion_age_test.R
# ==============================================================================
# Purpose:
#   Analyse the preliminary stem-lesion experiment comparing disease responses
#   in six- and nine-week-old quinoa plants following inoculation with
#   Sclerotinia sclerotiorum.
#
# Input:
#   data/phenotypic_data.xlsx
#   Worksheet: 2
#
# Outputs:
#   Descriptive statistics, Kruskal-Wallis tests, Dunn's post hoc comparisons,
#   and Wilcoxon tests for the stem-lesion age experiment
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

library(readxl)
library(dplyr)


# ==============================================================================
# 2. Define input file
# ==============================================================================

phenotypic_file <- file.path(
  "data",
  "phenotypic_data.xlsx"
)


# ==============================================================================
# 3. Import stem-lesion age-test data
# ==============================================================================

stem_lesion_data <- read_excel(
  phenotypic_file,
  sheet = 2
)


# ==============================================================================
# 4. Prepare data for analysis
# ==============================================================================

stem_test_data <- stem_lesion_data %>%
  transmute(
    Accession = factor(trimws(as.character(Accession))),
    Week = trimws(as.character(Week)),
    LesionSize = as.numeric(as.character(LesionSize))
  ) %>%
  filter(
    !is.na(Accession),
    !is.na(Week),
    !is.na(LesionSize),
    Week %in% c("6 Weeks", "9 Weeks")
  ) %>%
  mutate(
    Week = factor(
      Week,
      levels = c("6 Weeks", "9 Weeks")
    )
  ) %>%
  droplevels()

# Check factor levels
levels(stem_test_data$Accession)
levels(stem_test_data$Week)

# Check sample numbers
stem_test_data %>%
  count(Week, Accession)
str(stem_test_data)
# Check descriptive statistics
descriptive_results <- stem_test_data %>%
  group_by(Week, Accession) %>%
  summarise(
    n = n(),
    mean = mean(LesionSize),
    sd = sd(LesionSize),
    median = median(LesionSize),
    minimum = min(LesionSize),
    maximum = max(LesionSize),
    .groups = "drop"
  )

print(descriptive_results)

# Six-week plants
kw_6 <- kruskal.test(
  LesionSize ~ Accession,
  data = subset(stem_test_data, Week == "6 Weeks")
)

# Nine-week plants
kw_9 <- kruskal.test(
  LesionSize ~ Accession,
  data = subset(stem_test_data, Week == "9 Weeks")
)

kw_6
kw_9
kruskal_results <- data.frame(
  Week = c("6 Weeks", "9 Weeks"),
  Chi_square = c(
    unname(kw_6$statistic),
    unname(kw_9$statistic)
  ),
  df = c(
    unname(kw_6$parameter),
    unname(kw_9$parameter)
  ),
  p_value = c(
    kw_6$p.value,
    kw_9$p.value
  )
)

kruskal_results
dunn_6 <- FSA::dunnTest(
  LesionSize ~ Accession,
  data = subset(stem_test_data, Week == "6 Weeks"),
  method = "bh"
)$res

dunn_9 <- FSA::dunnTest(
  LesionSize ~ Accession,
  data = subset(stem_test_data, Week == "9 Weeks"),
  method = "bh"
)$res

dunn_6
dunn_9
significant_dunn_6 <- dunn_6 %>%
  filter(P.adj <= 0.05)

significant_dunn_9 <- dunn_9 %>%
  filter(P.adj <= 0.05)

significant_dunn_6
significant_dunn_9
wilcox_results <- stem_test_data %>%
  group_by(Accession) %>%
  group_modify(~ {
    
    values_6 <- .x$LesionSize[.x$Week == "6 Weeks"]
    values_9 <- .x$LesionSize[.x$Week == "9 Weeks"]
    
    test_result <- wilcox.test(
      x = values_9,
      y = values_6,
      paired = FALSE,
      exact = FALSE,
      alternative = "two.sided"
    )
    
    tibble(
      n_6_weeks = length(values_6),
      n_9_weeks = length(values_9),
      median_6_weeks = median(values_6),
      median_9_weeks = median(values_9),
      W = unname(test_result$statistic),
      p_value = test_result$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_adjusted_BH = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adjusted_BH <= 0.001 ~ "***",
      p_adjusted_BH <= 0.01  ~ "**",
      p_adjusted_BH <= 0.05  ~ "*",
      TRUE                   ~ "ns"
    )
  )

print(wilcox_results, n = Inf)
kruskal_results
dunn_6
dunn_9
wilcox_results


# ==============================================================================
# 5. Save results
# ==============================================================================

# Create workbook
wb <- createWorkbook()

#--------------------------------------------------
# Kruskal-Wallis
#--------------------------------------------------
addWorksheet(wb, "Kruskal_Wallis")
writeData(wb, "Kruskal_Wallis", kruskal_results)

#--------------------------------------------------
# Dunn test - 6 weeks
#--------------------------------------------------
addWorksheet(wb, "Dunn_6Weeks")
writeData(wb, "Dunn_6Weeks", dunn_6)

#--------------------------------------------------
# Dunn test - 9 weeks
#--------------------------------------------------
addWorksheet(wb, "Dunn_9Weeks")
writeData(wb, "Dunn_9Weeks", dunn_9)

#--------------------------------------------------
# Wilcoxon
#--------------------------------------------------
addWorksheet(wb, "Wilcoxon")

writeData(
  wb,
  "Wilcoxon",
  as.data.frame(wilcox_results)
)

#--------------------------------------------------
# Descriptive statistics (recommended)
#--------------------------------------------------
addWorksheet(wb, "Descriptive")

writeData(
  wb,
  "Descriptive",
  descriptive_results
)





# ==============================================================================
# 6. Plotting
# ==============================================================================

# Use the same accession order as the detached-leaf plot
accession_order <- c(
  "Bouchane-1",
  "Nde-09",
  "BO-11",
  "RU-5",
  "NL-6",
  "PUC-mix-red"
)

stem_lesion_data <- stem_lesion_data %>%
  dplyr::mutate(
    Accession = factor(
      as.character(Accession),
      levels = accession_order
    ),
    Week = factor(
      Week,
      levels = c("6 Weeks", "9 Weeks")
    )
  )

letters_df <- letters_df %>%
  dplyr::mutate(
    Accession = factor(
      as.character(Accession),
      levels = accession_order
    ),
    Week = factor(
      Week,
      levels = c("6 Weeks", "9 Weeks")
    )
  )


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p_stem_lesion <- ggplot(
  stem_lesion_data,
  aes(
    x = Accession,
    y = LesionSize
  )
) +
  
  # Individual observations
  geom_point(
    position = position_jitter(
      width = 0.10,
      height = 0,
      seed = 123
    ),
    shape = 21,
    size = 1.6,
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
  
  # Mean
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 2.6,
    colour = "black"
  ) +
  
  # Dunn-test letters
  geom_text(
    data = letters_df,
    aes(
      x = Accession,
      y = y_pos,
      label = Letter
    ),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "bold",
    colour = "black",
    size = 3
  ) +
  
  facet_wrap(
    ~Week,
    nrow = 1,
    scales = "fixed",
    drop = TRUE,
    labeller = as_labeller(
      c(
        "6 Weeks" = "6-week-old plants",
        "9 Weeks" = "9-week-old plants"
      )
    )
  ) +
  
  scale_x_discrete(
    limits = accession_order,
    drop = TRUE
  ) +
  
  scale_y_continuous(
    limits = c(0, NA),
    breaks = scales::breaks_pretty(n = 6),
    expand = expansion(
      mult = c(0, 0.10)
    )
  ) +
  
  labs(
    x = "Accession",
    y = "Lesion length (cm)"
  ) +
  
  theme_classic(
    base_size = 9,
    base_family = "Arial"
  ) +
  
  theme(
    legend.position = "none",
    
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
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      colour = "black",
      size = 9,
      face = "bold"
    ),
    
    panel.spacing = grid::unit(
      6,
      "mm"
    ),
    
    plot.margin = margin(
      4, 4, 4, 4,
      unit = "mm"
    )
  )

p_stem_lesion


ggsave(
  filename = paste0(
    "Stem_lesion_length_6_and_9_weeks.tiff"
  ),
  plot = p_stem_lesion,
  device = "tiff",
  width = 174,
  height = 110,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# ==============================================================================
#  END OF THE ANALYSIS
# ==============================================================================