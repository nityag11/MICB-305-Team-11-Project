#**Indicator Taxa Analysis**

# Step 1: Loading R Libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggpubr)
library(microbiome)
library(indicspecies)
library(writexl)
library(ANCOMBC)

# Step 2: Loading Data
metadata = read.csv('ms_metadata_corrected.csv', row.names = 1) 
taxonomy = read.delim('taxonomy.tsv', row.names = 1) 
counts = read.delim('feature-table.txt', skip = 1, row.names = 1, check.names = F)

# Step 3: Making a Phyloseq Object
taxonomy_formatted = taxonomy |>
  separate(col = Taxon, 
           into = c('Domain','Phylum','Class','Order','Family','Genus','Species'),
           sep=';', fill='right') |>
  select(-Confidence) |>
  as.matrix()

counts_formatted = counts |>
  as.matrix()

meta_names = read.csv('ms_metadata_corrected.csv', row.names = 1) |>
  names()

meta_data = read.csv('ms_metadata_corrected.csv', row.names = 1,
                     skip = 1,
                     header = F)

table(metadata[1,]== meta_data[1,])

names(meta_data) = meta_names

ps = phyloseq(sample_data(meta_data),
              otu_table(counts_formatted, taxa_are_rows = T),
              tax_table(taxonomy_formatted))

# Step 4: Running Indicator Taxa Analysis on "administration" Column 
# ---> [focus on "oral" vs. "infused"])
# ---> at Genus level

ps_genus = tax_glom(ps,'Genus')
ps_relab = microbiome::transform(ps_genus, 'compositional')
ps_filt = filter_taxa(ps_relab, function(x) mean(x) > 0.0001, TRUE)
otu_table = data.frame(otu_table(ps_filt))

set.seed(1)
indval = multipatt(t(otu_table), 
                    cluster= ps_filt@sam_data$administration,
                    control = how(nperm = 999))

summary(indval, indvalcomp = TRUE)


# ---> extracting it into a table
indval_table = as.data.frame(indval$sign)

genus_to_plot = indval_table |> 
  filter(stat > 0.7 | p.value < 0.05) |> 
  rownames()

genus_to_plot


df_of_taxa = prune_taxa(genus_to_plot, ps_filt) |> 
  psmelt()

df_of_taxa |> 
  ggplot(aes(administration, Abundance, fill = administration)) +
  geom_boxplot(position = position_dodge(width = 1.5),outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2, alpha = 0.15) +
  theme_classic(base_size = 10) +
  facet_wrap(~Phylum, ncol = 3, scales = 'free')

df_of_taxa
