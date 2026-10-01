# MANUSCRIPT ANALYSIS ====

# Belowground carbon stocks are shaped by edaphic heterogeneity
# and tree community structure in an eastern Amazonian rainforest
#
# Author: Anthony Barbosa
# Affiliation: Federal University of Pará, Brazil
# Contact: barbosaanthony.eco@gmail.com
#
# Last update: 2026-10-01


# "_______________________________________________________________"

# This script contains the statistical analyses and figures used in the
# manuscript. The workflow follows the original analysis script, with
# organization, documentation, and reproducibility improvements.

# 1. PACKAGES ==== 

# Install packages manually if they are not already available.
# Do NOT use install.packages() inside this analysis script.

library(lmerTest)
library(lme4)
library(DHARMa)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(vegan)
library(cowplot)
library(corrplot)
library(ggrepel)
library(scales)

# 2. DATA IMPORT ====
# In the original workflow:
#   - "ab" contains plot-level data.
#   - "plot" contains subplot-level data.
#
# The original script imported both datasets from the clipboard. This is kept
# here as an option for reproducing the analysis interactively in RStudio.

# Original import workflow:
# ab = read.table("clipboard", header = TRUE)
# plot = read.table("clipboard", header = TRUE)

# Example file-based import:
# ab = read.csv("data/plot_level_data.csv", header = TRUE)
# plot = read.csv("data/subplot_level_data.csv", header = TRUE)

# In the original files, ab columns are: PLOT | SITE | SOM | SMB | P | P_Rem | Fe | Al | pH | Clay | Richness | Diversity | Density | Basal_area

# SOM: Soil Organic Matter (g/kg)
# SMB: Soil Microbial Biomass (g/kg)
# P: Phosphorus (mg/kg)
# P_Rem: Remnant Phosphorus (mg/kg)
# Fe: Iron (mg/kg)
# Al: Aluminum (mg/kg)
# pH: Soil pH
# Clay: Clay content (%)
# Richness: Tree species richness
# Diversity: Tree species diversity (Shannon index)
# Density: Tree density (individuals per hectare)
# Basal_area: Tree basal area (m²/ha)

# At the subplot level, the plot columns are: PLOT | SUBPLOT | Richness | Density | Fine_roots | Coarse_roots | Litter | SITE

# 3. LINEAR MIXED-EFFECTS MODELS ====

# These models evaluate the relationships between:
#   - fine-root biomass,
#   - litter biomass, and
#   - coarse-root biomass
#
# and tree species richness (Riqueza) and site (Local).
#
# Plot is included as a random intercept to account for the hierarchical
# structure of the subplot-level data.
#
# The first two responses are modeled using linear mixed-effects models.
# Coarse-root biomass is modeled using a Gamma distribution with a log link.

fine_roots_raw = lmer(
  fine_roots ~ Richness + SITE + (1 | Plot),
  data = plot
)

summary(fine_roots_raw)


litter_raw = lmer(
  Litter ~ Richness + SITE + (1 | Plot),
  data = plot
)

summary(litter_raw)

coarse_roots = glmer(
  Coarse_roots ~ Richness + SITE + (1 | Plot),
  family = Gamma(link = "log"),
  data = plot
)

summary(coarse_roots)

# 4. MODEL DIAGNOSTICS ====

# DHARMa is used to assess model residuals and potential deviations from model
# assumptions. Diagnostics should be examined for each fitted model.

res_fine_roots = simulateResiduals(fine_roots_raw)
plot(res_fine_roots)

res_litter = simulateResiduals(litter_raw)
plot(res_litter)

res_coarse_roots = simulateResiduals(coarse_roots)
plot(res_coarse_roots)

# 5. MODEL SUMMARY AND EFFECT-SIZE PLOT ====

# broom.mixed is used to extract fixed-effect estimates and confidence
# intervals from the mixed-effects models.
library(broom.mixed)

#' [Extract fixed effects from each model]

# This function extracts the fixed effects, 95% confidence intervals, and
# p-values from each model.
#
# The intercept is excluded because the figure is intended to display the
# effects of the explanatory variables.

tidy_model = function(model, model_name) {
  
  tidy(model, effects = "fixed", conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      Model = model_name,
      Significance = ifelse(
        p.value < 0.05,
        "Significant",
        "Not significant"
      )
    )
}

# Combine results from the three models into a single data frame.
df_plot = bind_rows(
  tidy_model(d, "Fine roots"),
  tidy_model(e, "Litter"),
  tidy_model(c, "Coarse roots")
)


# Standardize variable labels for the figure

# The model output uses the original variable names from the dataset.
# Here, they are translated into English for presentation in the manuscript.

df_plot_en = df_plot %>%
  mutate(
    term = case_when(
      term == "Riqueza" ~ "Richness",
      term == "LocalSN" ~ "SiteSN",
      term == "LocalIGB" ~ "SiteIGB",
      TRUE ~ term
    )
  )

#' [Create effect-size plot]

# Points represent model estimates and horizontal lines represent 95%
# confidence intervals.

# The vertical dashed line at zero indicates no estimated effect.

grap = ggplot(
  df_plot_en,
  aes(
    x = estimate,
    y = term,
    xmin = conf.low,
    xmax = conf.high,
    color = Significance
  )
) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.5,
    linetype = "dashed",
    color = "gray30"
  ) +
  
  # Horizontal confidence intervals.
  # Using geom_errorbar() with orientation = "y" avoids relying on the
  # deprecated geom_errorbarh() function.
  geom_errorbar(
    orientation = "y",
    height = 0.22,
    linewidth = 1.2
  ) +
  
  # Model estimates.
  geom_point(
    aes(fill = Significance),
    shape = 21,
    size = 4.5,
    stroke = 1.4
  ) +
  
  # Each response variable is shown in a separate panel.
  facet_wrap(
    ~ Model,
    scales = "free_y"
  ) +
  
  # Significant estimates are filled, whereas non-significant estimates
  # remain white.
  scale_fill_manual(
    values = c(
      "Significant" = "#2C7FB8",
      "Not significant" = "white"
    )
  ) +
  
  scale_color_manual(
    values = c(
      "Significant" = "#2C7FB8",
      "Not significant" = "#D7301F"
    )
  ) +
  
  labs(
    x = "Estimated effect (β)",
    y = NULL,
    color = NULL,
    fill = NULL
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    text = element_text(
      family = "Arial",
      size = 12
    ),
    strip.text = element_text(
      family = "Arial",
      size = 12,
      face = "bold"
    ),
    axis.title.x = element_text(
      family = "Arial",
      size = 12,
      face = "bold"
    ),
    axis.text.x = element_text(
      family = "Arial",
      size = 12,
      face = "bold"
    ),
    axis.text.y = element_text(
      family = "Arial",
      size = 12,
      face = "bold"
    ),
    axis.title.y = element_text(
      family = "Arial",
      face = "bold"
    )
  )


# Display the figure in RStudio.
grap

# 6. PERMANOVA, PERMDISP, AND PCoA (Multivariate view) ====

# Convert plot and site identifiers to factors.
ab$PLOT = as.factor(ab$PLOT)
ab$SITE = as.factor(ab$SITE)

#' ´[Standardize plot-level environmental and vegetation variables]

# The variables in columns 3:17 (ab datafile, section 2) are standardized before calculating
# Euclidean distances.

# The standardized variables are organized into three groups:
#   - columns 3:10  = environmental variables
#   - columns 11:14 = vegetation structure
#   - columns 15:17 = organic input

ab1 = ab
ab1[3:17] = decostand(
  ab1[3:17],
  method = "standardize"
)

#' [Figure theme]

theme_paper = theme_classic(base_size = 8) +
  theme(
    text = element_text(
      family = "Arial",
      size = 12
    ),
    plot.title = element_text(
      face = "bold",
      size = 12
    ),
    axis.title = element_text(
      face = "bold",
      size = 12
    ),
    axis.text = element_text(
      size = 12
    ),
    legend.text = element_text(
      size = 12
    ),
    legend.title = element_text(
      size = 12
    )
  )


#' [Function for PERMANOVA, PERMDISP, and PCoA]

# This function:

# 1. Selects the variables specified by 'cols'.
# 2. Calculates a Euclidean distance matrix.
# 3. Tests differences among sites using PERMANOVA.
# 4. Tests differences in multivariate dispersion using PERMDISP.
# 5. Performs a Principal Coordinates Analysis (PCoA).
# 6. Calculates the centroid of each site.
# 7. Adds convex hulls around site groups.
# 8. Produces a PCoA plot displaying both PERMANOVA and PERMDISP results.

# The function returns the figure together with the statistical results,
# allowing the numerical results to be inspected separately in R.


run_permanova_plot = function(
    data,
    cols,
    group_var,
    title_plot,
    show_legend = TRUE
) {
  
  # Select variables and grouping factor.
  sub = data[, cols]
  group = data[[group_var]]
  
  # Calculate Euclidean distances among plots.
  dist_mat = vegdist(
    sub,
    method = "euclidean"
  )
  
  # PERMANOVA
  # --------------------------------------------------------------------------
  
  # Test whether multivariate composition differs between sites.
  #
  # 9,999 permutations are used to obtain the permutation-based p-value.
  
  perm <- adonis2(
    dist_mat ~ group,
    permutations = 9999
  )
  
  
  # --------------------------------------------------------------------------
  # PERMDISP
  # --------------------------------------------------------------------------
  
  # PERMDISP tests whether differences in multivariate dispersion may
  # contribute to the PERMANOVA result.
  
  disp <- betadisper(
    dist_mat,
    group
  )
  
  disp_anova <- anova(disp)
  
  
  # --------------------------------------------------------------------------
  # Principal Coordinates Analysis (PCoA)
  # --------------------------------------------------------------------------
  
  pcoa <- cmdscale(
    dist_mat,
    eig = TRUE,
    k = 2
  )
  
  ord <- as.data.frame(pcoa$points)
  
  colnames(ord) <- c(
    "PCoA1",
    "PCoA2"
  )
  
  ord$Group <- group
  
  
  # --------------------------------------------------------------------------
  # Group centroids
  # --------------------------------------------------------------------------
  
  centers <- ord %>%
    group_by(Group) %>%
    summarise(
      C1 = mean(PCoA1),
      C2 = mean(PCoA2),
      .groups = "drop"
    )
  
  
  # Add centroid coordinates to each observation.
  ord2 <- left_join(
    ord,
    centers,
    by = "Group"
  )
  
  
  # --------------------------------------------------------------------------
  # Convex hulls
  # --------------------------------------------------------------------------
  
  # Convex hulls are used only to visually delimit the distribution of
  # observations within each site.
  
  hull <- ord %>%
    group_by(Group) %>%
    slice(chull(PCoA1, PCoA2))
  
  
  # --------------------------------------------------------------------------
  # Statistical labels
  # --------------------------------------------------------------------------
  
  perm_label <- paste0(
    "PERMANOVA: Pseudo F = ",
    round(perm$F[1], 3),
    "; p = ",
    round(perm$`Pr(>F)`[1], 4)
  )
  
  
  disp_label <- paste0(
    "PERMDISP: F = ",
    round(disp_anova$`F value`[1], 3),
    "; p = ",
    round(disp_anova$`Pr(>F)`[1], 4)
  )
  
  
  # Site colors used consistently across the three panels.
  
  MyColor <- c(
    "maroon4",
    "royalblue"
  )
  
  
  # --------------------------------------------------------------------------
  # PCoA figure
  # --------------------------------------------------------------------------
  
  p <- ggplot(
    ord2,
    aes(
      PCoA1,
      PCoA2,
      color = Group
    )
  ) +
    
    # Convex hulls.
    geom_polygon(
      data = hull,
      aes(fill = Group),
      alpha = 0.15,
      color = NA
    ) +
    
    # Lines connecting observations to their group centroid.
    geom_segment(
      aes(
        xend = C1,
        yend = C2
      ),
      color = "gray60"
    ) +
    
    # Individual plots.
    geom_point(
      size = 3
    ) +
    
    # Group centroids.
    geom_point(
      data = centers,
      aes(
        C1,
        C2
      ),
      size = 5
    ) +
    
    scale_color_manual(
      values = MyColor
    ) +
    
    scale_fill_manual(
      values = MyColor
    ) +
    
    theme_paper +
    
    labs(
      title = title_plot,
      x = "PCoA1",
      y = "PCoA2"
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 19
      )
    ) +
    
    # Statistical results are displayed in the upper-right corner.
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = perm_label,
      hjust = 1,
      vjust = 2,
      size = 4.5
    ) +
    
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = disp_label,
      hjust = 1,
      vjust = 4,
      size = 4.5
    )
  
  
  # Remove the legend when requested.
  
  if (!show_legend) {
    p <- p +
      theme(
        legend.position = "none"
      )
  }
  
  
  # Return the plot and the statistical results.
  
  return(
    list(
      plot = p,
      permanova = perm,
      permdisp = disp_anova
    )
  )
}

#' [RUN MULTIVARIATE ANALYSES]

# Environmental variables: columns 3:10.
res_amb = run_permanova_plot(
  ab1,
  3:10,
  "SITE",
  "Environmental variables",
  show_legend = TRUE
)


# Vegetation structure: columns 11:14.
res_veg = run_permanova_plot(
  ab1,
  11:14,
  "SITE",
  "Vegetation structure",
  show_legend = TRUE
)


# Organic input: columns 15:17.
res_org = run_permanova_plot(
  ab1,
  15:17,
  "SITE",
  "Organic input",
  show_legend = TRUE
)

# Keep the legend in the organic-input panel so it can be extracted
# and displayed once in the combined figure.

res_org$plot = res_org$plot +
  theme(
    legend.position = "right"
  )

#' [COMBINE THE THREE PCoA PANELS]

legend = get_legend(
  res_org$plot
)

# Remove the legend from the individual panels.
p1 = res_amb$plot +
  theme(legend.position = "none")

p2 = res_veg$plot +
  theme(legend.position = "none")

p3 = res_org$plot +
  theme(legend.position = "none")

# Arrange the three panels horizontally.
plots_row = plot_grid(
  p1,
  p2,
  p3,
  ncol = 3,
  align = "hv",
  labels = c("A", "B", "C"),
  label_fontfamily = "Arial"
)

# Add the shared legend to the right.
panel = plot_grid(
  plots_row,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.12)
)

panel

# 7. UNIVARIATE COMPARISONS BETWEEN SITES ====

# The analysis compares each numeric variable between the two study sites.
#
# For each variable, normality is evaluated separately within each site.
# When the normality criterion is met in both groups, a two-sample t-test is
# performed. Otherwise, a Wilcoxon rank-sum test is used.

# Because multiple variables are tested, p-values are adjusted using the
# Benjamini-Hochberg (BH) procedure to control the false discovery rate (FDR).

#' [Identify numeric variables]-------------------------------------------------

vars = names(
  ab1
)[
  sapply(
    ab1,
    is.numeric
  )
]

#' [Run site comparisons]

# Create an empty data frame to store the results.
resultados_df = data.frame(
  Variable = vars,
  p_value = NA_real_,
  Test = NA_character_
)

for (i in seq_along(vars)) {
  
  v <- vars[i]
  
  
  # --------------------------------------------------------------------------
  # Test normality within each site
  # --------------------------------------------------------------------------
  
  # Shapiro-Wilk tests are performed separately for each site.
  
  normal <- by(
    ab1[[v]],
    ab1$SITE,
    function(x) shapiro.test(x)$p.value > 0.05
  )
  
  
  # --------------------------------------------------------------------------
  # Select the statistical test
  # --------------------------------------------------------------------------
  
  # If the normality criterion is met in both sites, use a two-sample t-test.
  # Otherwise, use the Wilcoxon rank-sum test.
  
  if (all(normal)) {
    
    test_result = t.test(
      as.formula(
        paste(v, "~ SITE")
      ),
      data = ab1
    )
    
    resultados_df$Test[i] = "t-test"
    
  } else {
    
    test_result = wilcox.test(
      as.formula(
        paste(v, "~ SITE")
      ),
      data = ab1
    )
    
    resultados_df$Test[i] = "Wilcoxon"
  }
  
  
# Store the raw p-value.
  resultados_df$p_value[i] = test_result$p.value
}
# Inspect the results before multiple-testing correction.
resultados_df


#' [Multiple-testing correction]

# The Benjamini-Hochberg procedure controls the expected false discovery rate
# when multiple statistical tests are performed.

resultados_df$p_adjust_BH = p.adjust(
  resultados_df$p_value,
  method = "BH"
)

# Add an indicator for statistically significant results after BH correction.
resultados_df$Significant_BH = resultados_df$p_adjust_BH < 0.05

# Display the final table.
resultados_df

# 8. PEARSON CORRELATIONS ====

# Calculate the Pearson correlation matrix for the environmental and biological
# variables included in columns 3 to 17.

cor_data = ab1[3:17]

cor_mat = cor(
  cor_data,
  method = "pearson"
)
cor_mat

#' [Correlation plot]

corrplot(
  cor_mat,
  method = "color",
  type = "upper",
  tl.cex = 0.8,
  addCoef.col = "black",
  number.cex = 0.7
)

# 9. PRINCIPAL COMPONENT ANALYSIS ====

# The PCA summarizes multivariate variation among plots based on the variables
# included in columns 3 to 17.
#
# Variables are standardized to zero mean and unit variance before the PCA so
# that variables measured on different scales contribute comparably.

#' [Select variables]
vars = ab[, 3:17]

#' [Run PCA]

# Standardization is performed directly within prcomp().
# This avoids scaling the data twice.
pca = prcomp(
  vars,
  center = TRUE,
  scale. = TRUE
)

summary(pca)

#' [PCA scores]

scores = as.data.frame(pca$x)
scores$SITE = ab$SITE[as.numeric(rownames(scores))]
scores$PLOT = ab$PLOT[as.numeric(rownames(scores))]

#' [PCA loadings]

loadings = as.data.frame(pca$rotation)
loadings$Variable = rownames(loadings)

#' [Explained variance

var_exp = (pca$sdev^2) / sum(pca$sdev^2)
var_pc1 = round(var_exp[1] * 100, 2)
var_pc2 = round(var_exp[2] * 100, 2)

# Contribution of variables to PC1 and PC2

loadings = loadings %>%
  mutate(
    contrib_PC1 = PC1^2 / sum(PC1^2) * 100,
    contrib_PC2 = PC2^2 / sum(PC2^2) * 100
  )

# Identify the three variables with the largest contributions to each axis

top_PC1 = loadings %>%
  arrange(desc(contrib_PC1)) %>%
  slice(1:3) %>%
  pull(Variable)

top_PC2 = loadings %>%
  arrange(desc(contrib_PC2)) %>%
  slice(1:3) %>%
  pull(Variable)

# Inspect the variables contributing most strongly to each axis.
top_PC1
top_PC2

#' [Classify variables for visualization]

# Variables occurring among the top three for both axes are assigned to the
# axis for which they have the larger contribution.

loadings = loadings %>%
  mutate(
    highlights = case_when(
      
      # Top three variables for PC1 only
      Variable %in% top_PC1 &
        !(Variable %in% top_PC2) ~ "PC1",
      
      # Top three variables for PC2 only
      Variable %in% top_PC2 &
        !(Variable %in% top_PC1) ~ "PC2",
      
      # Variables among the top three for both axes
      Variable %in% top_PC1 &
        Variable %in% top_PC2 &
        contrib_PC1 >= contrib_PC2 ~ "PC1",
      
      Variable %in% top_PC1 &
        Variable %in% top_PC2 &
        contrib_PC2 > contrib_PC1 ~ "PC2",
      
      # All remaining variables
      TRUE ~ "Others"
    )
  )

# Inspect variable contributions and classification.
loadings %>%
  select(
    Variable,
    contrib_PC1,
    contrib_PC2,
    Highlights
  ) %>%
  arrange(desc(contrib_PC1))

# Scale loading vectors for visualization
mult = min(
  
  (max(scores$PC1) - min(scores$PC1)) /
    (max(loadings$PC1) - min(loadings$PC1)),
  
  (max(scores$PC2) - min(scores$PC2)) /
    (max(loadings$PC2) - min(loadings$PC2))
  
) * 0.7

####################################################################