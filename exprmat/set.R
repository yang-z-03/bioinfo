
parser <- argparse::ArgumentParser(
  prog = "set",
  description = "set variable"
)

parser $ add_argument(
  "-f", type = "character", dest = "f", default = "",
  help = "the object data file, stored in .rds"
)

parser $ add_argument(
  "obj", type = "character",
  help = "the object name in the table"
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

shared[[pargs $ obj]] <- readRDS(pargs $ f)
