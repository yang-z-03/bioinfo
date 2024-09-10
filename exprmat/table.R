
parser <- argparse::ArgumentParser(
  prog = "table",
  description = "iterate all possible values"
)

parser $ add_argument(
  type = "character", dest = "name",
  help = "object to view"
)

parser $ add_argument(
  "-c", type = "character", dest = "col", default = "",
  help = "if the object is a dataframe, you may also specify column"
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

if (shared[[pargs $ name]] |> is.data.frame()) {
  shared[[pargs $ name]][[pargs $ col]] |> table() |> print()
} else {
  shared[[pargs $ name]] |> table() |> print()
}
