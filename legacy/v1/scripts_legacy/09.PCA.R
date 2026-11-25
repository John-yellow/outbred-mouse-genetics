#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
if (is.na(script_file) || !nzchar(script_file)) {
  project_root <- normalizePath(".", mustWork = FALSE)
} else {
  project_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
}

get_arg <- function(flag, default) {
  hit <- which(args == flag)
  if (length(hit) > 0 && hit[1] < length(args)) return(args[hit[1] + 1])
  default
}

resolve_path <- function(path) {
  if (grepl("^/", path)) path else file.path(project_root, path)
}

pca_prefix <- resolve_path(get_arg("--prefix", Sys.getenv("PCA_PREFIX", "results/genetics/obm_pca")))
output_dir <- resolve_path(get_arg("--output-dir", Sys.getenv("FIGURE_DIR", "results/figures")))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(ggplot2))

eigenvec_file <- paste0(pca_prefix, ".eigenvec")
eigenval_file <- paste0(pca_prefix, ".eigenval")
if (!file.exists(eigenvec_file)) stop("PCA eigenvector file not found: ", eigenvec_file)
if (!file.exists(eigenval_file)) stop("PCA eigenvalue file not found: ", eigenval_file)

eigenvec <- read.table(eigenvec_file, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
if (ncol(eigenvec) < 4) stop("PCA eigenvector file must contain FID, IID, and at least two PCs")
colnames(eigenvec) <- c("group", "IID", paste0("PC", seq_len(ncol(eigenvec) - 2)))

normalise_group <- function(x) {
  x <- as.character(x)
  x[grepl("^(ICR1|ICR_H)", x)] <- "ICR1"
  x[grepl("^(ICR2|ICR_W)", x)] <- "ICR2"
  x[grepl("^(KM|K)", x)] <- "KM"
  x[x == "I"] <- "ICR1"
  x
}
eigenvec$group <- normalise_group(eigenvec$group)

eigenval <- scan(eigenval_file, quiet = TRUE)
if (length(eigenval) < 2) stop("PCA eigenvalue file contains fewer than two values")
explained_variance <- round(eigenval / sum(eigenval) * 100, 2)

colors <- c(ICR1 = "coral", ICR2 = "#D3B356", KM = "cadetblue3")
PCAplot <- ggplot(eigenvec, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = group), size = 3, alpha = 0.8) +
  labs(
    title = NULL,
    x = paste0("PC1 (", explained_variance[1], "% variance)"),
    y = paste0("PC2 (", explained_variance[2], "% variance)")
  ) +
  scale_color_manual(values = colors, na.value = "grey50") +
  theme_minimal() +
  theme(
    text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  )

ggsave(file.path(output_dir, "PCA.pdf"), plot = PCAplot, width = 7, height = 6.5)
save.image(file = file.path(output_dir, "PCA.RData"))
