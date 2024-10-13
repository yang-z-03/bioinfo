
parser <- argparse::ArgumentParser(
  prog = "normp",
  description = "normalize quantitative proteomic data"
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}
