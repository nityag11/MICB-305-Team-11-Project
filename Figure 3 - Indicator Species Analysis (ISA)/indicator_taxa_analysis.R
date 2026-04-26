title: "indicator_species_analysis"
author: "juliana"
date: "2026-04-25"
---

# Indicator Species Analysis

#Indicator Species Analysis (ISA) was used to identify bacterial taxa dominant or specific in treatment groups with or without probiotic supplementation.

## Loading R Libraries

```{r}
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggpubr)
library(microbiome)
library(indicspecies)
library(writexl)
```

## Loading and Filtering Dataset

```{r}
#Loading in Dataset Files
taxonomy = read.delim("taxonomy.tsv", row.names = 1, check.names = FALSE)
tree = read_tree("tree.nwk")
counts = read.delim("feature-table.txt", skip = 1, row.names = 1, check.names = FALSE)
metadata = read.csv("ms_metadata_curated.csv", row.names = 1)


#Filtering Metadata
metadata_filt = metadata |> 
  filter(probiotics != "NA")

metadata_filt$group <- ifelse(metadata_filt$treatments_grouped == "Untreated", "Untreated", paste(metadata_filt$treatments_grouped, metadata_filt$probiotics, sep = "_"))

metadata_omit <- metadata_filt |> 
  filter(group != "Untreated", group != "Control_No", group != "Control_Yes")


#Creating the Phyloseq Object
taxonomy_formatted = taxonomy |>
  separate(col = Taxon,
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep =";", fill = "right") |>
  select(-Confidence) |>
  as.matrix()

counts_formatted = counts |>
  as.matrix()

ps = phyloseq(sample_data(metadata_omit), 
              otu_table(counts_formatted, taxa_are_rows = T),
              tax_table(taxonomy_formatted),
              tree)

```
## Running ISA on Filtered Metadata

```{r}
ps_genus = tax_glom(ps,'Genus') #ASVs aggregated to Genus level
ps_relab = microbiome::transform(ps_genus, 'compositional') #Phyloseq converted to relative abundance
ps_filt = filter_taxa(ps_relab, function(x) mean(x) > 0.0001, TRUE) #Abundance filter applied to keep taxa whose mean relative abundance is greater than 0.01%
otu_table = data.frame(otu_table(ps_filt))

set.seed(1)
indval = multipatt(t(otu_table), 
                    cluster= ps_filt@sam_data$group, #Samples were clustered by treatment group
                    control = how(nperm = 999))

#Summary Table of Results
summary_table = summary(indval, indvalcomp = TRUE)
```
## Extracting Results into a Table

```{r}
#Index number values were converted into their respective group names
indval_table = as.data.frame(indval$sign) |> 
  mutate(index = case_when(
    index == 1 ~ "Non-Oral (-)", 
    index == 2 ~ "Non-Oral (+)", 
    index == 3 ~ "Oral (-)", 
    index == 4 ~ "Oral (+)", 
    index == 5 ~ "Non-Oral Group (-/+)", 
    index == 6 ~ "No probiotics", 
    index == 7 ~ "Not applicable", #Number not present in table
    index == 8 ~ "Not applicable", #Number not present in table
    index == 9 ~ "Non-Oral (+) & Oral (+)", 
    index == 10 ~ "Oral Group (-/+)", 
    index == 11 ~ "All except Oral (+)", 
    index == 12 ~ "All except Oral (-)", 
    index == 13 ~ "All except Non-Oral (+)", 
    index == 14 ~ "All except Non-Oral (-)", 
    index == 15 ~ "All groups"
  ))
  

genus_to_plot = indval_table |> 
  filter(p.value <= 0.05, stat >= 0.5) |> #Significance: p-value less than 0.05, stat value greater than 0.5
  rownames()

indval_join = indval_table |> 
  rownames_to_column("OTU")

df_of_taxa = prune_taxa(genus_to_plot, ps_filt) |> 
  psmelt() |> 
  left_join(indval_join) |> 
  mutate(admin_probiotics = case_when(group == "Non-Oral_No" ~ "Non-Oral (-)", 
                                    group == "Non-Oral_Yes" ~ "Non-Oral (+)", 
                                    group == "Oral_Yes" ~ "Oral (+)",
                                    group == "Oral_No" ~ "Oral (-)")) |> #Renaming group names for consistency 
  group_by(admin_probiotics, Genus) |>
  mutate(mean_abundance = mean(Abundance)) |>
  ungroup()

```
## Graphing Results on a Dot Plot

#Note: The final graph used in Figure 3 was manually edited using Microsoft Paint for clarity and aesthetics. These edits include: (1) Renaming the x-axis "admin_probiotics" to "Probiotics". (2) Renaming "index" in the legend to "Indicator Treatment Group(s)". (3) Renaming "mean_abundance" in the legend to "Mean Abundance". (4) Adding "Non-Oral" and "Oral" grouping labels to the top of the graph. (5) Adding "Administration" label to the top of the graph. (6) Replacing the labels below the graph ("Non-Oral (-), Non-Oral (+), Oral (-), Oral (+)") with "(-)" and "(+)" to match new labels at the top of the graph.

```{r}
indicator_plot = df_of_taxa |> 
  ggplot(aes(x = admin_probiotics, y = Genus)) +
  geom_point(aes(colour = index, size = mean_abundance), alpha = 0.5) +
  theme(element_text(size=16), 
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank())

indicator_plot
