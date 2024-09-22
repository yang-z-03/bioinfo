
parser <- argparse::ArgumentParser(
  prog = "intgmeta",
  description = "integrate gene metadata"
)

parser $ add_argument(
  "-i", dest = "input", type = "character", nargs = "*", default = c(),
  help = paste("input sample directory name")
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

if (length(pargs $ input) > 0) {
  run_names <- pargs $ input
} else {
  run_names <- list.dirs(".", full.names = FALSE, recursive = FALSE)
}

merged_gene_info <- NULL
for (runname in run_names) {
  if (runname == "data") next
  if (runname == "norm") next
  
  geneinfo <- readRDS(
    paste(".", runname, "norm", "genes-meta.rds", sep = "/")
  )

  if (merged_gene_info |> is.null()) {
    merged_gene_info <- geneinfo
  } else {
    merged_gene_info <- dplyr::bind_rows(merged_gene_info, geneinfo)
  }
}

dup <- duplicated(pull(merged_gene_info, "seurat_names")) |
  pull(merged_gene_info, "seurat_names") == ""

merged_gene_info <- merged_gene_info[!dup, ]

# filtering and ordering

ord <- c()
for (cx in rownames(shared[["seurat"]])) {
  id <- which(merged_gene_info $ seurat_names == cx)
  ord <- c(ord, id[1])
}

meta <- merged_gene_info[ord, ]
saveRDS(meta, "norm/genes-meta.rds")
shared[["meta_gene"]] <- meta
