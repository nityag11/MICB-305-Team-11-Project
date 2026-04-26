suppressPackageStartupMessages({
  library(phyloseq)
  library(pheatmap)
  library(stringr)
})


indicator_genera <- c(
  "Anaerotruncus",
  "Actinomyces",
  "Anaerosporobacter",
  "[Eubacterium]_fissicatena_group",
  "Frisingicoccus",
  "Bacillus",
  "Staphylococcus",
  "Dielma",
  "Flavonifractor",
  "Pediococcus",
  "Anaerococcus"
)

min_mean_relative_abundance <- 1e-4
target_genus_count <- 20
heatmap_file <- "heatmap.png"
matrix_file <- "correlation_matrix.tsv"
selected_genera_file <- "selected_genera.tsv"
run_summary_file <- "aim2_5_run_summary.tsv"

read_feature_table <- function(path) {
  all_lines <- readLines(path, warn = FALSE)
  data_lines <- all_lines[!startsWith(all_lines, "# Constructed from biom file")]

  if (length(data_lines) < 2) {
    stop("feature-table.txt does not contain enough data after removing comment lines.")
  }

  data_lines[1] <- sub("^#OTU ID", "Feature.ID", data_lines[1])
  data_lines[1] <- sub("^\t", "Feature.ID\t", data_lines[1])

  table_df <- read.delim(
    text = paste(data_lines, collapse = "\n"),
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE
  )

  rownames(table_df) <- table_df[[1]]
  table_df[[1]] <- NULL
  table_df
}

extract_rank <- function(taxon_string, prefix) {
  rank_match <- str_match(taxon_string, paste0(prefix, "([^;]+)"))
  rank_values <- rank_match[, 2]
  rank_values <- str_trim(rank_values)
  rank_values[rank_values == ""] <- NA
  rank_values
}

find_script_dir <- function() {
  script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_args) > 0) {
    script_path <- sub("^--file=", "", script_args[1])
    return(dirname(normalizePath(script_path)))
  }

  candidates <- c(
    getwd(),
    file.path(getwd(), "aim2_5_analysis"),
    "D:/Wiz LXBB/LEO-0417-03/aim2_5_analysis"
  )

  for (candidate in candidates) {
    candidate <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
    if (file.exists(file.path(candidate, "305.R"))) {
      return(candidate)
    }
  }

  stop("Could not determine the script directory. Set the working directory to aim2_5_analysis and rerun.")
}

script_dir <- find_script_dir()
project_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
setwd(script_dir)

feature <- read_feature_table(file.path(project_dir, "feature-table.txt"))

taxonomy <- read.delim(
  file.path(project_dir, "taxonomy.tsv"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

metadata <- read.delim(
  file.path(project_dir, "ms_metadata_corrected.tsv"),
  header = TRUE,
  sep = "\t",
  fill = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)

feature_id_col <- if ("Feature.ID" %in% names(taxonomy)) "Feature.ID" else "Feature ID"
sample_id_col <- if ("sample_name" %in% names(metadata)) "sample_name" else "sample-id"

taxonomy$Kingdom <- extract_rank(taxonomy$Taxon, "d__")
taxonomy$Phylum <- extract_rank(taxonomy$Taxon, "p__")
taxonomy$Class <- extract_rank(taxonomy$Taxon, "c__")
taxonomy$Order <- extract_rank(taxonomy$Taxon, "o__")
taxonomy$Family <- extract_rank(taxonomy$Taxon, "f__")
taxonomy$Genus <- extract_rank(taxonomy$Taxon, "g__")
taxonomy$Species <- extract_rank(taxonomy$Taxon, "s__")

taxonomy_index <- match(rownames(feature), taxonomy[[feature_id_col]])
matched_taxonomy <- taxonomy[taxonomy_index, , drop = FALSE]
keep_rows <- !is.na(matched_taxonomy$Genus)

feature <- feature[keep_rows, , drop = FALSE]
matched_taxonomy <- matched_taxonomy[keep_rows, , drop = FALSE]

common_samples <- intersect(colnames(feature), metadata[[sample_id_col]])
if (length(common_samples) < 3) {
  stop("Fewer than 3 overlapping samples were found between the feature table and metadata.")
}

feature <- feature[, common_samples, drop = FALSE]
metadata <- metadata[match(common_samples, metadata[[sample_id_col]]), , drop = FALSE]
rownames(metadata) <- metadata[[sample_id_col]]

taxonomy_matrix <- cbind(
  Kingdom = matched_taxonomy$Kingdom,
  Phylum = matched_taxonomy$Phylum,
  Class = matched_taxonomy$Class,
  Order = matched_taxonomy$Order,
  Family = matched_taxonomy$Family,
  Genus = paste0("g__", matched_taxonomy$Genus),
  Species = matched_taxonomy$Species
)
rownames(taxonomy_matrix) <- rownames(feature)

ps <- phyloseq(
  otu_table(as.matrix(feature), taxa_are_rows = TRUE),
  sample_data(metadata),
  tax_table(as.matrix(taxonomy_matrix))
)

ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = TRUE)
otu_genus <- as(otu_table(ps_genus), "matrix")
if (!taxa_are_rows(ps_genus)) {
  otu_genus <- t(otu_genus)
}

sample_totals <- colSums(otu_genus)
keep_samples <- sample_totals > 0
otu_genus <- otu_genus[, keep_samples, drop = FALSE]
sample_totals <- sample_totals[keep_samples]

if (ncol(otu_genus) < 3) {
  stop("Fewer than 3 samples remain after removing zero-count samples.")
}

genus_names <- as.character(tax_table(ps_genus)[, "Genus"])
genus_names <- sub("^g__", "", genus_names)
rownames(otu_genus) <- genus_names

genus_rel_abundance <- sweep(otu_genus, 2, sample_totals, "/")
genus_mean_abundance <- rowMeans(genus_rel_abundance)
genus_variance <- apply(genus_rel_abundance, 1, var)

candidate_genera <- rownames(genus_rel_abundance)[
  genus_mean_abundance >= min_mean_relative_abundance & genus_variance > 0
]

priority_genera <- intersect(indicator_genera, candidate_genera)
remaining_genera <- setdiff(candidate_genera, priority_genera)
remaining_genera <- remaining_genera[order(genus_mean_abundance[remaining_genera], decreasing = TRUE)]
selected_genera <- unique(c(priority_genera, remaining_genera))[seq_len(min(target_genus_count, length(unique(c(priority_genera, remaining_genera)))))]

if (length(selected_genera) < 2) {
  stop("Fewer than 2 indicator genera remained after filtering.")
}

selected_rel_abundance <- genus_rel_abundance[selected_genera, , drop = FALSE]
selected_rel_abundance <- selected_rel_abundance[
  apply(selected_rel_abundance, 1, var) > 0,
  ,
  drop = FALSE
]

if (nrow(selected_rel_abundance) < 2) {
  stop("Fewer than 2 indicator genera with non-zero variance are available for correlation analysis.")
}

correlation_matrix <- cor(
  t(selected_rel_abundance),
  method = "spearman",
  use = "pairwise.complete.obs"
)

display_genera <- paste0("g__", rownames(selected_rel_abundance))
rownames(correlation_matrix) <- display_genera
colnames(correlation_matrix) <- display_genera

write.table(
  data.frame(Genus = rownames(selected_rel_abundance)),
  file = selected_genera_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    metric = c(
      "samples_initial_overlap",
      "samples_removed_zero_counts",
      "samples_final",
      "candidate_genera",
      "selected_genera"
    ),
    value = c(
      length(common_samples),
      sum(!keep_samples),
      ncol(selected_rel_abundance),
      length(candidate_genera),
      nrow(selected_rel_abundance)
    )
  ),
  file = run_summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  correlation_matrix,
  file = matrix_file,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

heatmap_plot <- pheatmap(
  correlation_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45
)

png(heatmap_file, width = 3000, height = 2400, res = 300)
grid::grid.draw(heatmap_plot$gtable)
dev.off()
