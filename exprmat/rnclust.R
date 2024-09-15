
parser <- argparse::ArgumentParser(
  prog = "rnclust",
  description = "rename (or merge) the cluster values"
)

parser $ add_argument(
  "-c", dest = "cluster", type = "character", default = "",
  help = paste("specify the cluster identification")
)

parser $ add_argument(
  "-n", dest = "names",
  type = "character", default = c(), nargs = "+",
  help = paste("names to assign each cluster")
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

# rename cluster columns, including those generated with sc3 and seurat.
# this can be viewed with lsclust and defclust.

cls <- shared[["seurat"]][[pargs $ cluster]] |> class()

if (cls[1] == "factor") {
  labels(shared[["seurat"]][[pargs $ cluster]]) <- pargs $ names
  labels(pull(shared[["meta_sample"]], pargs $ cluster)) <- pargs $ names
}
