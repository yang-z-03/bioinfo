
parser <- argparse::ArgumentParser(
  prog = "lsclust",
  description = "list the cluster values"
)

parser $ add_argument(
  "-c", dest = "cluster", type = "character", default = "",
  help = paste("specify the cluster identification")
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

cls <- shared[["seurat"]][[pargs $ cluster]] |> class()

if (cls[1] == "factor") {
  print(shared[["seurat"]][[pargs $ cluster]] |> table())
}
