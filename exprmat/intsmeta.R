
parser <- argparse::ArgumentParser(
  prog = "intsmeta",
  description = "integrate sample metadata"
)

parser $ add_argument(
  "-c", dest = "col", type = "character", nargs = "*", default = c(),
  help = paste("column names you would like to extract")
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

meta <- data.frame(
  id = colnames(shared[["seurat"]])
)

for (x in pargs $ col) {
  meta[[x]] <- shared[["seurat"]] @ meta.data |> pull(x)
}

meta <- tibble(meta)
saveRDS(meta, "norm/samples-meta.rds")
shared[["meta_sample"]] <- meta
