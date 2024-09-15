
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

cls <- shared[["seurat"]] @ meta.data |> pull(pargs $ cluster) |> class()

if (cls[1] == "factor") {
  cat(blue("active ident:"))
  print(Idents(shared[["seurat"]]) |> table())
  cat(crlf)
  cat(blue("selected ident:"))
  print(shared[["seurat"]] @ meta.data |> pull(pargs $ cluster) |> table())
}
