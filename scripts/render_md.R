library(getopt)

script_dir <- if (exists("snakemake")) {
  dirname(snakemake@script)
} else if (rstudioapi::isAvailable()) {
  dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  dirname(normalizePath(commandArgs(trailingOnly = FALSE) |>
                          grep("--file=", x = _, value = TRUE) |>
                          sub("--file=", "", x = _)))
}

spec <- matrix(c(
  "datafolder", "d", 1, "character",
  "samples", "s", 1, "character",
  "help", "h", 0, "logical"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

template <- paste0(script_dir, "/rqc.Rmd")

folder       <- normalizePath(opt$datafolder, mustWork = TRUE)

if(!is.null(opt$samples)) {
  samples_file <- normalizePath(opt$samples, mustWork = TRUE)
} else {
  samples_file <- paste0(folder, "/reference/samples.tsv")
}

samples_name <- gsub(".tsv", "", basename(samples_file))

if(!dir.exists(paste0(folder, "/report"))) dir.create(paste0(folder, "/report"), recursive = TRUE)

rmarkdown::render(
      input = template,
      params = list(project = folder, samples_file = samples_file),
      output_file = paste0(folder, "/report/", samples_name,".report.html"))

