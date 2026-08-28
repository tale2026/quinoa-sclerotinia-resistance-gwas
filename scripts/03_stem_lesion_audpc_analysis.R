# ==============================================================================
# Script: 03_stem_lesion_audpc_analysis.R
# ==============================================================================
# Purpose:
#   Analyse lesion-length progression and area under the disease progress curve
#   (AUDPC) among quinoa accessions following stem inoculation with
#   Sclerotinia sclerotiorum.
#
# Input:
#   data/phenotypic_data.xlsx
#   Worksheet: 3
#
# Outputs:
#   Descriptive statistics, accession comparisons, correlation analyses,
#   and figures for lesion length and AUDPC
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
library(openxlsx)
library(GGally)
library(ggplot2)



# ==============================================================================
# 2. Define input file and import data
# ==============================================================================

phenotypic_file <- file.path(
  "data",
  "phenotypic_data.xlsx"
)


stem_lesion_raw <- read_excel(
  phenotypic_file,
  sheet = 3
)


# ==============================================================================
# 4. Prepare data for analysis
# ==============================================================================

stem_lesion_data <- stem_lesion_raw %>%
  transmute(
    accession = factor(trimws(as.character(accession))),
    id = factor(trimws(as.character(id))),
    across(
      c(day7, day12, day17, audpc),
      ~ as.numeric(as.character(.x))
    ),
    country = factor(trimws(as.character(country))),
    region = factor(trimws(as.character(region)))
  ) %>%
  filter(
    !is.na(accession),
    !is.na(id)
  ) %>%
  droplevels()

#inspect data
str(stem_lesion_data)
summary(stem_lesion_data)


## ============================================================================
## 5. Accession-wise descriptive statistics
## ============================================================================

desc_by_accession <- stem_lesion_data %>%
  group_by(accession) %>%
  summarise(
    across(
      c(day7, day12, day17, audpc),
      list(
        N      = ~sum(!is.na(.)),
        Mean   = ~mean(., na.rm = TRUE),
        SD     = ~sd(., na.rm = TRUE),
        Median = ~median(., na.rm = TRUE),
        Min    = ~min(., na.rm = TRUE),
        Max    = ~max(., na.rm = TRUE),
        Q1     = ~quantile(., 0.25, na.rm = TRUE),
        Q3     = ~quantile(., 0.75, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# Export descriptive statistics
write.xlsx(
  desc_by_accession,
  file = "descriptive_statistics_by_accession.xlsx",
  overwrite = TRUE
)

## ============================================================================
## 6. Calculate accession means
## ============================================================================

accession_means <- stem_lesion_data %>%
  group_by(accession) %>%
  summarise(
    day7  = mean(day7, na.rm = TRUE),
    day12 = mean(day12, na.rm = TRUE),
    day17 = mean(day17, na.rm = TRUE),
    audpc = mean(audpc, na.rm = TRUE),
    .groups = "drop"
  )

## ============================================================================
## 7. Pearson correlation analysis
## ============================================================================

cor_matrix <- cor(
  accession_means[, c("day7", "day12", "day17", "audpc")],
  method = "pearson"
)

print(cor_matrix)

## ============================================================================
## 8. Correlation plot
## ============================================================================

plot_data <- accession_means[, c("day7", "day12", "day17", "audpc")]

colnames(plot_data) <- c(
  "Lesion length\n(7 dpi)",
  "Lesion length\n(12 dpi)",
  "Lesion length\n(17 dpi)",
  "AUDPC"
)

cor_plot <- ggpairs(
  plot_data,
  upper = list(
    continuous = wrap(
      "cor",
      size = 5,
      prefix = "r = "
    )
  ),
  lower = list(
    continuous = wrap(
      "smooth",
      method = "lm",
      se = FALSE,
      alpha = 0.75,
      size = 0.4
    )
  ),
  diag = list(
    continuous = wrap(
      "densityDiag",
      alpha = 0.5
    )
  ),
  progress = FALSE
) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid = element_blank()
  )

cor_plot

## ============================================================================
## 9. Assumption testing
## 9.1 Normality of residuals (Shapiro-Wilk)
## ============================================================================

traits <- c("day7", "day12", "day17", "audpc")

for (trait in traits) {
  
  cat("\n=====================================\n")
  cat("Trait:", trait, "\n")
  cat("=====================================\n")
  
  # Fit linear model
  model <- lm(as.formula(paste(trait, "~ accession")), data = stem_lesion_data)
  
  # Shapiro-Wilk test on residuals
  print(shapiro.test(residuals(model)))
  
  # Histogram of residuals
  hist(
    residuals(model),
    main = paste("Residuals -", trait),
    xlab = "Residuals",
    col = "lightgrey",
    border = "black"
  )
  
  # QQ plot
  qqnorm(
    residuals(model),
    main = paste("Q-Q plot -", trait)
  )
  qqline(
    residuals(model),
    col = "red",
    lwd = 2
  )
}

library(FSA)

## ============================================================================
## 9.2 Homogeneity of variance (Levene's test)
## ============================================================================

library(car)

traits <- c("day7", "day12", "day17", "audpc")

for (trait in traits) {
  
  cat("\n=====================================\n")
  cat("Trait:", trait, "\n")
  cat("=====================================\n")
  
  formula <- as.formula(paste(trait, "~ accession"))
  
  print(leveneTest(formula, data = stem_lesion_data))
}

## ============================================================================
## 10. One-way ANOVA, Kruskal-Wallis and pairwise comparisons
## ============================================================================

library(openxlsx)
library(FSA)

traits <- c("day7", "day12", "day17", "audpc")

anova_results <- data.frame()
kw_results <- data.frame()

tukey_list <- list()
dunn_list  <- list()

for(trait in traits){
  
  formula <- as.formula(paste(trait, "~ accession"))
  
  ##----------------------------------------------------
  ## One-way ANOVA
  ##----------------------------------------------------
  
  model <- aov(formula, data = stem_lesion_data)
  
  aov_tab <- summary(model)[[1]]
  
  anova_results <- rbind(
    anova_results,
    data.frame(
      Trait = trait,
      DF_between = aov_tab[1,"Df"],
      DF_within  = aov_tab[2,"Df"],
      F_value    = round(aov_tab[1,"F value"],3),
      P_value    = aov_tab[1,"Pr(>F)"]
    )
  )
  
  ##----------------------------------------------------
  ## Tukey HSD
  ##----------------------------------------------------
  
  tukey <- TukeyHSD(model)
  
  tukey_df <- data.frame(
    Comparison = rownames(tukey$accession),
    tukey$accession,
    Trait = trait,
    row.names = NULL
  )
  
  tukey_list[[trait]] <- tukey_df
  
  ##----------------------------------------------------
  ## Kruskal-Wallis
  ##----------------------------------------------------
 str(stem_lesion_data)
  kw <- kruskal.test(formula, data = stem_lesion_data)
  
  kw_results <- rbind(
    kw_results,
    data.frame(
      Trait = trait,
      Chi_square = round(as.numeric(kw$statistic),3),
      DF = as.numeric(kw$parameter),
      P_value = kw$p.value
    )
  )
  
  ##----------------------------------------------------
  ## Dunn's test (Holm adjustment)
  ##----------------------------------------------------
  
  dunn <- dunnTest(
    formula,
    data = stem_lesion_data,
    method = "holm"
  )
  
  dunn_df <- dunn$res
  
  dunn_df$Trait <- trait
  
  dunn_list[[trait]] <- dunn_df
  
}

##----------------------------------------------------
## Export to Excel
##----------------------------------------------------

wb <- createWorkbook()

addWorksheet(wb, "ANOVA")
writeData(wb, "ANOVA", anova_results)

addWorksheet(wb, "Kruskal_Wallis")
writeData(wb, "Kruskal_Wallis", kw_results)

for(trait in traits){
  
  addWorksheet(wb, paste0("Tukey_", trait))
  writeData(
    wb,
    paste0("Tukey_", trait),
    tukey_list[[trait]]
  )
  
  addWorksheet(wb, paste0("Dunn_", trait))
  writeData(
    wb,
    paste0("Dunn_", trait),
    dunn_list[[trait]]
  )
  
}

saveWorkbook(
  wb,
  "ANOVA_KruskalWallis_Pairwise.xlsx",
  overwrite = TRUE
)

anova_results
kw_results

str(sheet3)


## ============================================================================
## 11. Highland vs. Lowland comparison (accession means)
## ============================================================================


accession_means$region <- factor(accession_means$region)

# Traits to analyse
traits <- c("day7", "day12", "day17", "audpc")

ttest_results <- data.frame()
wilcox_results <- data.frame()

for(trait in traits){
  
  formula <- as.formula(paste(trait, "~ region"))
  
  ## Welch's t-test
  tt <- t.test(formula, data = accession_means)
  
  ttest_results <- rbind(
    ttest_results,
    data.frame(
      Trait = trait,
      Mean_Highland = mean(accession_means[[trait]][accession_means$region=="highland"], na.rm=TRUE),
      Mean_Lowland  = mean(accession_means[[trait]][accession_means$region=="lowland"], na.rm=TRUE),
      t_value = round(as.numeric(tt$statistic),3),
      DF = round(as.numeric(tt$parameter),2),
      P_value = tt$p.value
    )
  )
  
  ## Wilcoxon rank-sum test
  wt <- wilcox.test(formula,
                    data = accession_means,
                    exact = FALSE)
  
  wilcox_results <- rbind(
    wilcox_results,
    data.frame(
      Trait = trait,
      W = as.numeric(wt$statistic),
      P_value = wt$p.value
    )
  )
}

## ============================================================================
## Export results
## ============================================================================

write.xlsx(
  list(
    Welch_t_test = ttest_results,
    Wilcoxon = wilcox_results
  ),
  file = "Region_comparison.xlsx",
  overwrite = TRUE
)

ttest_results
wilcox_results


## ============================================================================
## AUDPC by country (accession means)
## ============================================================================
library(dplyr)
library(ggplot2)

# Summary statistics
audpc_summary <- stem_lesion_data %>%
  group_by(accession, country) %>%
  summarise(
    mean = mean(audpc, na.rm = TRUE),
    sd   = sd(audpc, na.rm = TRUE),
    .groups = "drop"
  )

audpc_summary <- stem_lesion_data %>%
  group_by(accession, country) %>%
  summarise(
    mean = mean(audpc, na.rm = TRUE),
    sd = sd(audpc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(country, accession)

accession_levels <- audpc_summary$accession

stem_lesion_data$accession <- factor(
  stem_lesion_data$accession,
  levels = accession_levels
)

audpc_summary$accession <- factor(
  audpc_summary$accession,
  levels = accession_levels
)

#changed the order of accessions

accession_order <- stem_lesion_data %>%
  group_by(country, accession) %>%
  summarise(mean_audpc = mean(audpc, na.rm = TRUE), .groups = "drop") %>%
  arrange(country, mean_audpc) %>%  # lowest to highest within country
  pull(accession) %>%
  as.character()

stem_lesion_data <- stem_lesion_data %>%
  mutate(accession = factor(accession, levels = accession_order))

audpc_summary <- audpc_summary %>%
  mutate(accession = factor(accession, levels = accession_order))

aggregate(audpc ~ region, data = stem_lesion_data, median)

#------------------------------------------------------------
# 12. Broad-sense heritability (H²) for AUDPC
#------------------------------------------------------------

str(stem_lesion_data)
library(lme4)

# Fit random effects model
H2_model <- lmer(audpc ~ (1 | accession), data = stem_lesion_data)

# Extract variance components
vc <- as.data.frame(VarCorr(H2_model))

genetic_variance  <- vc$vcov[vc$grp == "accession"]
residual_variance <- vc$vcov[vc$grp == "Residual"]

# Average number of biological replicates per accession
n_rep <- mean(table(stem_lesion_data$accession))

# Entry-mean broad-sense heritability
H2 <- genetic_variance /
  (genetic_variance + (residual_variance / n_rep))

# Print results
cat("Genetic variance :", genetic_variance, "\n")
cat("Residual variance:", residual_variance, "\n")
cat("Average replicates:", round(n_rep, 2), "\n")
cat("Broad-sense heritability (H²):", round(H2, 4), "\n")


# Load necessary library
if (!require("lme4")) install.packages("lme4")
library(lme4)

# 1. Fit the mixed-effects model
# We treat 'accession' as a random effect to estimate its variance component.
# If 'region' or 'country' are experimental design factors (blocks), 
# include them as fixed or random effects.
model <- lmer(audpc ~ (1 | accession), data = stem_lesion_data)

# 2. Extract variance components
var_comp <- as.data.frame(VarCorr(model))
var_g <- var_comp[var_comp$grp == "accession", "vcov"]  # Genetic variance
var_e <- var_comp[var_comp$grp == "Residual", "vcov"]   # Residual (environmental) variance

# 3. Calculate Heritability (Entry-mean basis)
# Determine the average number of replicates (n) per accession
n_reps <- mean(table(stem_lesion_data$accession))

h2_entry_mean <- var_g / (var_g + (var_e / n_reps))

# 4. Calculate Heritability (Individual-plant basis)
h2_individual <- var_g / (var_g + var_e)

# Results
cat("Broad-sense Heritability (Entry-mean):", round(h2_entry_mean, 3), "\n")
cat("Broad-sense Heritability (Individual):", round(h2_individual, 3), "\n")



#------------------------------------------------------------
# 13.Plotting
#------------------------------------------------------------

# ------------------------------------------------------------
#  Function to extract the first valid metadata value
# ------------------------------------------------------------

first_valid_value <- function(x) {
  
  x <- trimws(as.character(x))
  
  x[
    x == "" |
      tolower(x) %in% c("na", "n/a", "nan")
  ] <- NA_character_
  
  valid_values <- x[!is.na(x)]
  
  if (length(valid_values) == 0) {
    return(NA_character_)
  }
  
  valid_values[1]
}


# ------------------------------------------------------------
# Extract country and region per accession
# ------------------------------------------------------------

accession_metadata <- stem_lesion_data %>%
  dplyr::mutate(
    accession = trimws(as.character(accession))
  ) %>%
  dplyr::filter(
    accession != "CHEN-425",
    !is.na(accession),
    accession != ""
  ) %>%
  dplyr::group_by(accession) %>%
  dplyr::summarise(
    country = first_valid_value(country),
    region  = first_valid_value(region),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    country = dplyr::if_else(
      is.na(country) | trimws(country) == "",
      "Unknown",
      trimws(as.character(country))
    ),
    
    region = dplyr::if_else(
      is.na(region) | trimws(region) == "",
      "Not specified",
      trimws(as.character(region))
    )
  )


# ------------------------------------------------------------
# Prepare raw AUDPC data
# ------------------------------------------------------------

audpc_plot_data <- stem_lesion_data %>%
  dplyr::transmute(
    accession = trimws(as.character(accession)),
    audpc = as.numeric(audpc)
  ) %>%
  dplyr::filter(
    accession != "CHEN-425",
    !is.na(accession),
    accession != "",
    !is.na(audpc)
  ) %>%
  dplyr::left_join(
    accession_metadata,
    by = "accession"
  ) %>%
  dplyr::mutate(
    region = dplyr::case_when(
      grepl("high", region, ignore.case = TRUE) ~ "Highland",
      grepl("low", region, ignore.case = TRUE)  ~ "Lowland",
      TRUE                                      ~ "Not specified"
    )
  )


# Confirm removal of CHEN-425
stopifnot(!"CHEN-425" %in% audpc_plot_data$accession)


# ------------------------------------------------------------
# Calculate accession mean and standard deviation
# ------------------------------------------------------------

audpc_summary_plot <- audpc_plot_data %>%
  dplyr::group_by(
    country,
    accession
  ) %>%
  dplyr::summarise(
    n = sum(!is.na(audpc)),
    mean_audpc = mean(audpc, na.rm = TRUE),
    sd_audpc = stats::sd(audpc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    sd_audpc = dplyr::if_else(
      is.na(sd_audpc),
      0,
      sd_audpc
    )
  )


# ------------------------------------------------------------
# Set country order
# Unknown is placed last
# ------------------------------------------------------------

known_countries <- audpc_summary_plot %>%
  dplyr::filter(country != "Unknown") %>%
  dplyr::distinct(country) %>%
  dplyr::arrange(country) %>%
  dplyr::pull(country)

country_order <- c(
  known_countries,
  "Unknown"
)


# ------------------------------------------------------------
# Order accessions by country and mean AUDPC
# ------------------------------------------------------------

accession_positions <- audpc_summary_plot %>%
  dplyr::mutate(
    country_order_number = match(
      country,
      country_order
    )
  ) %>%
  dplyr::arrange(
    country_order_number,
    mean_audpc,
    accession
  ) %>%
  dplyr::mutate(
    x_pos = dplyr::row_number()
  ) %>%
  dplyr::select(
    country,
    accession,
    x_pos
  )


# Add numerical x positions to the raw data
audpc_plot_data <- audpc_plot_data %>%
  dplyr::left_join(
    accession_positions,
    by = c("country", "accession")
  )


# Add numerical x positions to the summary data
audpc_summary_plot <- audpc_summary_plot %>%
  dplyr::left_join(
    accession_positions,
    by = c("country", "accession")
  )


# ------------------------------------------------------------
# Calculate positions of country brackets
# ------------------------------------------------------------

data_upper_limit <- max(
  audpc_plot_data$audpc,
  audpc_summary_plot$mean_audpc +
    audpc_summary_plot$sd_audpc,
  na.rm = TRUE
)

# Round upward to the nearest multiple of five
data_upper_limit <- ceiling(data_upper_limit / 5) * 5

bracket_spacing <- max(
  2,
  data_upper_limit * 0.055
)

tick_length <- max(
  0.8,
  data_upper_limit * 0.025
)


country_brackets <- accession_positions %>%
  dplyr::group_by(country) %>%
  dplyr::summarise(
    xmin = min(x_pos) - 0.38,
    xmax = max(x_pos) + 0.38,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    country_order_number = match(
      country,
      country_order
    )
  ) %>%
  dplyr::arrange(country_order_number) %>%
  dplyr::mutate(
    country_number = dplyr::row_number(),
    
    # Two alternating bracket levels
    bracket_y =
      data_upper_limit +
      bracket_spacing +
      ((country_number - 1) %% 2) * bracket_spacing,
    
    country_label = dplyr::case_when(
      country == "United Kingdom" ~ "UK",
      country == "United States"  ~ "USA",
      TRUE                        ~ as.character(country)
    )
  )


# Move all country brackets and names downward
country_brackets <- country_brackets %>%
  dplyr::mutate(
    bracket_y = bracket_y - 4
  )

country_brackets <- country_brackets %>%
  dplyr::mutate(
    label_angle = dplyr::if_else(
      country_label %in% c(
        "Bolivia",
        "Chile",
        "Peru",
        "USA",
        "Unknown"
      ),
      0,
      90
    ),
    
    label_hjust = dplyr::if_else(
      label_angle == 0,
      0.5,
      0
    ),
    
    label_vjust = dplyr::if_else(
      label_angle == 0,
      0,
      0.5
    )
  )
# Upper plot limit, including space for vertical country names
plot_upper_limit <- max(
  country_brackets$bracket_y,
  na.rm = TRUE
) + bracket_spacing * 2.8


# ------------------------------------------------------------
# X-axis labels
# ------------------------------------------------------------

x_axis_breaks <- accession_positions$x_pos
x_axis_labels <- accession_positions$accession


# ------------------------------------------------------------
# Create plot
# ------------------------------------------------------------

audpc_country_plot <- ggplot2::ggplot() +
  
  # Individual biological replicates
  ggplot2::geom_point(
    data = audpc_plot_data,
    ggplot2::aes(
      x = x_pos,
      y = audpc,
      shape = region
    ),
    position = ggplot2::position_jitter(
      width = 0.13,
      height = 0,
      seed = 123
    ),
    size = 0.95,
    stroke = 0.18,
    colour = "grey45",
    alpha = 0.75
  ) +
  
  # Mean ± standard deviation
  ggplot2::geom_errorbar(
    data = audpc_summary_plot,
    ggplot2::aes(
      x = x_pos,
      ymin = pmax(0, mean_audpc - sd_audpc),
      ymax = mean_audpc + sd_audpc
    ),
    width = 0.12,
    linewidth = 0.25,
    colour = "black"
  ) +
  
  # Mean
  ggplot2::geom_point(
    data = audpc_summary_plot,
    ggplot2::aes(
      x = x_pos,
      y = mean_audpc
    ),
    shape = 18,
    size = 1.30,
    colour = "black"
  ) +
  
  # Horizontal country brackets
  ggplot2::geom_segment(
    data = country_brackets,
    ggplot2::aes(
      x = xmin,
      xend = xmax,
      y = bracket_y,
      yend = bracket_y
    ),
    linewidth = 0.25,
    colour = "black"
  ) +
  
  # Left side of country brackets
  ggplot2::geom_segment(
    data = country_brackets,
    ggplot2::aes(
      x = xmin,
      xend = xmin,
      y = bracket_y,
      yend = bracket_y - tick_length
    ),
    linewidth = 0.25,
    colour = "black"
  ) +
  
  # Right side of country brackets
  ggplot2::geom_segment(
    data = country_brackets,
    ggplot2::aes(
      x = xmax,
      xend = xmax,
      y = bracket_y,
      yend = bracket_y - tick_length
    ),
    linewidth = 0.25,
    colour = "black"
  ) +
  
  # Vertical country names immediately above brackets
  ggplot2::geom_text(
    data = country_brackets,
    ggplot2::aes(
      x = (xmin + xmax) / 2,
      y = bracket_y + 0.35,
      label = country_label,
      angle = label_angle,
      hjust = label_hjust,
      vjust = label_vjust
    ),
    family = "Arial",
    fontface = "bold",
    size = 1.8,
    colour = "black"
  ) +
  
  ggplot2::scale_shape_manual(
    values = c(
      "Highland"      = 16,
      "Lowland"       = 17,
      "Not specified" = 15
    ),
    name = "Region"
  ) +
  
  ggplot2::scale_x_continuous(
    breaks = x_axis_breaks,
    labels = x_axis_labels,
    limits = c(
      0.3,
      max(accession_positions$x_pos) + 0.7
    ),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  
  ggplot2::scale_y_continuous(
    limits = c(0, plot_upper_limit),
    breaks = seq(
      0,
      45,
      by = 5
    ),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  
  ggplot2::labs(
    x = "Accession",
    y = "AUDPC"
  ) +
  
  ggplot2::guides(
    shape = ggplot2::guide_legend(
      title.position = "left",
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        size = 2,
        colour = "grey35",
        alpha = 1
      )
    )
  ) +
  
  ggplot2::theme_bw(
    base_size = 9,
    base_family = "Arial"
  ) +
  
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(
      size = 9,
      margin = ggplot2::margin(t = 6)
    ),
    
    axis.title.y = ggplot2::element_text(
      size = 9,
      margin = ggplot2::margin(r = 5)
    ),
    
    axis.text.x = ggplot2::element_text(
      size = 3.8,
      angle = 60,
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
    
    # Remove complete box
    panel.border = ggplot2::element_blank(),
    
    # Retain only left and bottom axes
    axis.line.x = ggplot2::element_line(
      colour = "black",
      linewidth = 0.50
    ),
    
    axis.line.y = ggplot2::element_line(
      colour = "black",
      linewidth = 0.50
    ),
    
    panel.grid.major.x = ggplot2::element_blank(),
    
    panel.grid.minor = ggplot2::element_blank(),
    
    panel.grid.major.y = ggplot2::element_line(
      colour = "grey85",
      linewidth = 0.35
    ),
    
    legend.position = "bottom",
    
    legend.title = ggplot2::element_text(
      size = 8,
      face = "bold"
    ),
    
    legend.text = ggplot2::element_text(
      size = 8
    ),
    
    plot.margin = ggplot2::margin(
      t = 4,
      r = 3,
      b = 3,
      l = 3,
      unit = "mm"
    )
  )


# ------------------------------------------------------------
# Display plot
# ------------------------------------------------------------

audpc_country_plot


# ------------------------------------------------------------
# Save publication-quality files
# ------------------------------------------------------------

# TIFF: 600 dpi
ggplot2::ggsave(
  filename = file.path(
    output_directory,
    "AUDPC_one_panel_country_brackets6.tiff"
  ),
  plot = audpc_country_plot,
  device = "tiff",
  width = 174,
  height = 120,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

------------------------------------------------------------
# END OF THE ANALYSIS
# ------------------------------------------------------------