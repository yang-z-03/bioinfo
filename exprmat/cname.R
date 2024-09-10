
parser <- argparse::ArgumentParser(
  prog = "cname",
  description = "view column names"
)

parser $ add_argument(type = "character", dest = "object",
                      help = "variable object name to view")

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

colnames(shared[[pargs $ object]]) |> print()
