# Load in libraries

library(tidyverse)
library(phyloseq)

# Load Datasets
taxonomy = read.delim("exports/taxonomy.tsv", row.names = 1, check.names = FALSE)
tree = read_tree("exports/tree.nwk")
counts = read.delim("exports/feature-table.txt", skip = 1, row.names = 1, check.names = FALSE)
metadata = read.delim("ms_metadata_curated.tsv", row.names = 1)

# Wrangle Tables

# Taxonomy 

taxonomy_formatted = taxonomy %>%
  separate(col = Taxon,
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep =";", fill = "right") %>%
  select(-Confidence) %>%
  as.matrix()

# Counts

counts_formatted = counts %>% as.matrix()

# Metadata 

