
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
