
parser <- argparse::ArgumentParser(
  prog = "chassay",
  description = "change seurat assay"
)

parser $ add_argument(
  type = "character", dest = "assay",
  help = "seurat assay name"
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

if (shared[["is_norm"]]) {
  if (pargs $ assay %in% SeuratObject::Assays(shared[["seurat"]])) {
    Seurat::DefaultAssay(shared[["seurat"]]) <- pargs $ assay
  } else cat(red("no such assay in seurat object."))
} else cat("you need to normalize first.")
