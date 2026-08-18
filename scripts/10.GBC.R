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

q_prefix <- resolve_path(get_arg("--q-prefix", Sys.getenv("Q_PREFIX", "results/genetics/obm_pruned_dataset.3")))
fam_file <- resolve_path(get_arg("--fam", Sys.getenv("FAM_FILE", "results/genetics/obm_pruned_dataset.fam")))
output_dir <- resolve_path(get_arg("--output-dir", Sys.getenv("FIGURE_DIR", "results/figures")))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(pheatmap)
})

q_file <- if (grepl("\\.Q$", q_prefix)) q_prefix else paste0(q_prefix, ".Q")
if (!file.exists(q_file)) stop("ADMIXTURE Q file not found: ", q_file)
if (!file.exists(fam_file)) stop("PLINK FAM file not found: ", fam_file)

Q <- read.table(q_file, header = FALSE, stringsAsFactors = FALSE)
fam <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
if (nrow(Q) != nrow(fam)) stop("Q matrix and FAM file have different numbers of samples")
if (ncol(Q) < 3) stop("The Q matrix must contain at least three ancestry components")

k <- min(3, ncol(Q))
Q <- Q[, seq_len(k), drop = FALSE]
colnames(Q) <- c("KM1", "ICR1", "ICR2")[seq_len(k)]
Q$group <- as.character(fam$V1)
Q$group[Q$group == "KM"] <- "KM1"
Q$group[grepl("^ICR_H", Q$group)] <- "ICR1"
Q$group[grepl("^ICR_W", Q$group)] <- "ICR2"

target_order <- c("ICR1", "ICR2", "KM1")
Q <- Q %>%
  mutate(group = factor(group, levels = target_order)) %>%
  arrange(group)
unique_groups <- unique(Q$group)
group_sizes <- table(factor(Q$group, levels = unique_groups))

ICR1c <- adjustcolor("coral", alpha.f = 0.7)
ICR2c <- adjustcolor("#D3B356", alpha.f = 0.7)
KM1c <- adjustcolor("cadetblue3", alpha.f = 0.7)

pdf(file.path(output_dir, "ADMIXTURE.pdf"), width = 15, height = 6)
bar_positions <- barplot(
  t(as.matrix(Q[, seq_len(k), drop = FALSE])),
  col = c(KM1c, ICR1c, ICR2c)[seq_len(k)],
  border = NA,
  space = rep(0.1, nrow(Q)),
  axes = FALSE
)
axis(side = 2, cex.axis = 1.4, font = 2, lwd = 2)
group_ends <- cumsum(group_sizes)
group_starts <- c(1, head(group_ends, -1) + 1)
x_centers <- sapply(seq_along(group_sizes), function(i) {
  mean(bar_positions[group_starts[i]:group_ends[i]])
})
text(x = x_centers, y = -0.05, labels = unique_groups, srt = 0, adj = 0.5, cex = 1.6, font = 2, xpd = TRUE)
dev.off()

Q_export <- Q
Q_export$group <- as.character(Q_export$group)
for (column in names(Q_export)) {
  if (is.numeric(Q_export[[column]])) Q_export[[column]] <- sprintf("%.5f", Q_export[[column]])
}
write.table(Q_export, file.path(output_dir, "ADMIXTURE_sorted_Q.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

Q_group <- Q %>%
  group_by(group) %>%
  summarise(across(all_of(colnames(Q)[seq_len(k)]), mean), .groups = "drop") %>%
  mutate(group = as.character(group))
for (column in names(Q_group)) {
  if (is.numeric(Q_group[[column]])) Q_group[[column]] <- sprintf("%.5f", Q_group[[column]])
}
write.table(Q_group, file.path(output_dir, "GBC.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
