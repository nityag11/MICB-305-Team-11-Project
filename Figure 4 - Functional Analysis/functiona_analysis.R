# Manuscript Functional Analysis

```{r}
# install.packages("devtools")
# devtools::install_github("biobakery/maaslin2")
# devtools::install_github("cafferychen777/ggpicrust2")
# install.packages('MicrobiomeStat')
# BiocManager::install("KEGGREST")
# install.packages("GGally")
# install.packages("BiocManager")
# BiocManager::install("ANCOMBC")

library(tidyverse)
library(dplyr)
library(phyloseq)
library(ggpicrust2)
library(Maaslin2)
library(ggplot2)
```

## General Data Wrangling

```{r}
# General Data Wrangling
metadata_ms = readRDS('../Michael_Cen_Assignment_5/phyloseq_ms_taxonomy.rds') %>% .@sam_data %>% data.frame() |>
  rownames_to_column('sample_name')

# Re-labeling different groups into oral and non-oral
ms_treatment_group = metadata_ms |>
  mutate(treatment_group = case_when(
    administration == "Subcutaneous" ~ "Non-oral",
    administration == "Infused" ~ "Non-oral",
    administration == "Untreated" ~ "Non-oral",
    administration == "Control" ~ "Control",
    administration == "Oral" ~ "Oral",
  )) |>
  relocate(treatment_group, .after = administration)

# Filter out control
oralvsnonoral = ms_treatment_group |>
  filter(treatment_group %in% c("Oral", "Non-oral"))

# Remove NA
no_na_oralvsnonoral = oralvsnonoral |>
  filter(!is.na(probiotics))

# Creating two tabs for later purposes
pro_oral_ms_treatment_group = no_na_oralvsnonoral |>
  mutate(
    probiotics_label = ifelse(probiotics == 1, "probiotic", "nonprobiotic"),
    admin_label = ifelse(administration == "Oral", "oral", "nonoral")
  )

metacyc = read.delim("path_abun_unstrat.tsv", row.names = 1, check.names = FALSE)

# Ensuring same sample size between MetaCyc results and Metadata
common_samples_filt_prooral <- intersect(colnames(metacyc), pro_oral_ms_treatment_group$sample_name)
pro_oral_ms_treatment_group <- pro_oral_ms_treatment_group |>
  filter(sample_name %in% common_samples_filt_prooral)
metacyc_pro_oral_filt <- metacyc |>
  select(all_of(common_samples_filt_prooral))

# Setting the default background for future MaAslin2 Tests:
metadata_maaslin <- pro_oral_ms_treatment_group |>
  mutate(probiotics = factor(ifelse(probiotics == 1, "yes", "no"),
                             levels = c("no", "yes")),
         administration = factor(admin_label,
                                 levels = c("nonoral", "oral"))
  )
metacyc_t <- as.data.frame(t(metacyc_pro_oral_filt))

shared_samples <- intersect(rownames(metacyc_t), metadata_maaslin$sample_name)

metacyc_t <- metacyc_t[shared_samples, , drop = FALSE]
metadata_maaslin_filtered <- metadata_maaslin |>
  filter(sample_name %in% shared_samples)

```

## MaAslin2 Data Analysis code:

```{r}
# MaAslin2 analysis on correlation between administration route and Pathway abundant With respect to Probiotic Usage

# Probiotic Group
metadata_probiotic <- metadata_maaslin_filtered |>
  filter(probiotics == "yes") |>
  column_to_rownames("sample_name")

metacyc_probiotic <- metacyc_t[rownames(metadata_probiotic), , drop = FALSE]

maaslin_probiotic <- Maaslin2(
  input_data = metacyc_probiotic,
  input_metadata = metadata_probiotic,
  output = "maaslin_probiotic_only",
  fixed_effects = c("administration"),
  max_significance = 0.05
)

results_probiotic_all <- read_tsv("maaslin_probiotic_only/all_results.tsv")

# Non-probiotic
metadata_nonprobiotic <- metadata_maaslin_filtered |>
  filter(probiotics == "no") |>
  column_to_rownames("sample_name")

metacyc_nonprobiotic <- metacyc_t[rownames(metadata_nonprobiotic), , drop = FALSE]

maaslin_nonprobiotic <- Maaslin2(
  input_data = metacyc_nonprobiotic,
  input_metadata = metadata_nonprobiotic,
  output = "maaslin_nonprobiotic_only",
  fixed_effects = c("administration"),
  max_significance = 0.05
)

results_nonprobiotic_all <- read_tsv("maaslin_nonprobiotic_only/all_results.tsv")

# Interaction Term
metadata_interaction <- metadata_maaslin_filtered |>
  mutate(
    interaction = as.numeric(administration == "oral") *
                  as.numeric(probiotics == "yes")
  ) |>
  column_to_rownames("sample_name")

interaction_maaslin <- Maaslin2(
  input_data = metacyc_t,
  input_metadata = metadata_interaction,
  output = "maaslin_interaction_results",
  fixed_effects = c("administration", "probiotics", "interaction")
)

results_interact_all <- read.delim("maaslin_interaction_results/all_results.tsv")

results_interaction_only <- results_interact_all |>
  filter(metadata == "interaction") |>
  arrange(qval)

# Identifying sample size in each group
table(metadata_probiotic$administration)
table(metadata_nonprobiotic$administration)
```

## Volcano Plot Generation

```{r}
# Volcano Plot Generation

## Probiotic Volcano Plot
probiotic_maaslin_volc <- results_probiotic_all |>
  mutate(significance = if_else(qval < 0.05, "Significant", "NS")) |>
  ggplot(aes(x = coef, y = -log10(qval))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  labs(
    x = "coef (Log2FC)",
    y = "-log10(q-value)"
  ) +
       xlim(-2.5, 2.5) +
  ylim(0, 2) +
  theme_minimal(base_size = 20)
probiotic_maaslin_volc

## Non-probiotic Volcano Plot
 nonprobiotic_maaslin_volc <- results_nonprobiotic_all |>
  mutate(significance = if_else(qval < 0.05, "Significant", "NS")) |>
  ggplot(aes(x = coef, y = -log10(qval))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  labs(
    x = "coef (Log2FC)",
    y = "-log10(q-value)"
  ) +
     xlim(-1, 1) +
  ylim(0, 2) +
  theme_minimal(base_size = 20)
nonprobiotic_maaslin_volc

## Interaction Term Volcano Plot
results_interaction_only <- results_interact_all |>
  filter(metadata == "interaction") |>
  arrange(qval)
  
interaction_volc <- results_interaction_only |>
  mutate(sig = qval < 0.05) |>
  ggplot(aes(x = coef, y = -log10(qval))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  labs(
    x = "coef (Log2FC)",
    y = "-log10(q-value)"
  ) +
  xlim(-1, 1) +
  ylim(0, 2) +
  theme_minimal(base_size = 20)
interaction_volc

# Saving Graphs:
# ggsave("probiotic_volcano_plot.png", plot = probiotic_maaslin_volc)
# ggsave("nonprobiotic_volcano_plot.png", plot = nonprobiotic_maaslin_volc)
# ggsave("interaction_volcano_plot.png", plot = interaction_volc)
```

## PCA Graph Generation

```{r}
# Generating PCA Graph for Probiotic and Non-probiotic data subsets

## Non-probiotic Samples
nonprobiotic_pca <- pathway_pca(
  abundance = metacyc_pro_oral_filt,
  metadata = metadata_nonprobiotic,
  group = "administration"
)
nonprobiotic_pca

## Probiotic Samples
probiotic_pca <- pathway_pca(
  abundance = metacyc_pro_oral_filt,
  metadata = metadata_probiotic,
  group = "administration"
)
probiotic_pca

# Saving Graphs
# ggsave("nonprobiotic_pca.png", plot = nonprobiotic_pca)
# ggsave("probiotic_pca.png", plot = probiotic_pca)
```
