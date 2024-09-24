
parser <- argparse::ArgumentParser(
  prog = "cd",
  description = "change working directory"
)

parser $ add_argument(
  type = "character", dest = "dir",
  help = "working directory (absolute or relative)"
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

setwd(pargs $ dir)
cat("set working directory to", green(pargs $ dir) |> italic(), crlf)

shared[["is_reference_assigned"]] <- FALSE
shared[["is_loaded"]] <- FALSE
shared[["is_qc"]] <- FALSE
shared[["is_norm"]] <- FALSE

if (file.exists("genome.rds")) {
  shared[["is_reference_assigned"]] <- TRUE
}

# autoloads

if (file.exists("norm/genes-meta.rds") &&
    file.exists("norm/samples-meta.rds") &&
    file.exists("norm/seurat.rds")) {

  shared[["is_ready"]] <- TRUE
  shared[["meta_sample"]] <- readRDS("norm/samples-meta.rds")
  shared[["meta_gene"]] <- readRDS("norm/genes-meta.rds")
  shared[["seurat"]] <- readRDS("norm/seurat.rds")
  cat(yellow("all normalization preprocesses are ready."), crlf)
  stop()
} else {
  shared[["is_ready"]] <- FALSE
  shared[["meta_sample"]] <- NULL
  shared[["meta_gene"]] <- NULL
  shared[["seurat"]] <- NULL
}

if (file.exists("norm/linear.rds") &&
      file.exists("norm/log.rds") &&
      file.exists("norm/seurat.rds")) {

  shared[["is_norm"]] <- TRUE
  shared[["meta_sample"]] <- readRDS("norm/samples-meta.rds")
  shared[["meta_gene"]] <- readRDS("norm/genes-meta.rds")
  shared[["counts"]] <- readRDS("qc/matrix.rds")
  shared[["seurat"]] <- readRDS("norm/seurat.rds")
  cat(blue("autoload from qc/matrix.rds and norm/*"), crlf)
  cat(yellow("all normalization preprocesses are ready."), crlf)
  stop()

} else {
  shared[["is_norm"]] <- FALSE
  shared[["meta_sample"]] <- NULL
  shared[["meta_gene"]] <- NULL
  shared[["seurat"]] <- NULL
}

if (file.exists("qc/matrix.rds") &&
      file.exists("qc/genes-meta.rds") &&
      file.exists("qc/samples-meta.rds")) {

  shared[["is_qc"]] <- TRUE
  shared[["counts"]] <- readRDS("qc/matrix.rds")
  shared[["meta_sample_raw"]] <- readRDS("qc/samples-meta.rds")
  shared[["meta_gene_raw"]] <- readRDS("qc/genes-meta.rds")
  cat(blue("autoload from qc/*"), crlf)
  stop()

} else {
  shared[["is_qc"]] <- FALSE
  shared[["meta_sample_raw"]] <- NULL
  shared[["meta_gene_raw"]] <- NULL
  shared[["counts"]] <- NULL
}

if (file.exists("features/matrix.rds") &&
      file.exists("features/genes-meta.rds") &&
      file.exists("features/samples-meta.rds")) {

  shared[["is_loaded"]] <- TRUE
  shared[["counts"]] <- readRDS("features/matrix.rds")
  shared[["meta_sample_raw"]] <- readRDS("features/samples-meta.rds")
  shared[["meta_gene_raw"]] <- readRDS("features/genes-meta.rds")
  cat(blue("autoload from features/*"), crlf)
  stop()

} else {
  shared[["is_loaded"]] <- FALSE
  shared[["meta_sample_raw"]] <- NULL
  shared[["meta_gene_raw"]] <- NULL
  shared[["counts"]] <- NULL
}
