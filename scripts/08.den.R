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

input_dir <- resolve_path(get_arg("--input-dir", Sys.getenv("MAP_DIR", "results/genetics/populations")))
output_dir <- resolve_path(get_arg("--output-dir", Sys.getenv("FIGURE_DIR", "results/figures")))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(CMplot)
  library(dplyr)
})

read_map <- function(population) {
  map_file <- file.path(input_dir, paste0(population, ".map"))
  if (!file.exists(map_file)) stop("Map file not found: ", map_file)
  map <- fread(map_file, header = FALSE)
  if (ncol(map) < 4) stop("PLINK map file has fewer than four columns: ", map_file)
  map <- map %>% select(SNP = 2, Chromosome = 1, Position = 4)
  map$Chromosome <- as.character(map$Chromosome)
  map$Chromosome[map$Chromosome == "23"] <- "X"
  map
}

maps <- setNames(lapply(c("ICR1", "ICR2", "KM"), read_map), c("ICR1", "ICR2", "KM"))

old_wd <- getwd()
setwd(output_dir)

for (population in names(maps)) {
  CMplot(
    maps[[population]],
    plot.type = "d",
    bin.size = 1e6,
    chr.den.col = c("#6F9D87FF", "#FFD966", "#FF7878FF"),
    bin.breaks = seq(0, 15000, 3000),
    file = "pdf",
    file.output = TRUE,
    verbose = TRUE,
    file.name = population
  )
}
setwd(old_wd)
