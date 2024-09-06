
parser_qc <- argparse::ArgumentParser(prog = "qc")

parser_qc $ add_argument("-d", type = "integer", dest = "mindet",
                         help = "minimal genes detected per cell",
                         default = 500)

parser_qc $ add_argument("-c", type = "double", dest = "mincell",
                         help = "minimal expressing cells per gene (in percentage)", # nolint
                         default = 0.1)

parser_qc $ add_argument(
  "-m", type = "integer", dest = "mindet",
  help = paste("the upper limit of mitochondrial gene percentage. set to 0",
               "to automatically pick the threshold by isOutlier()"),
  default = 0
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser_qc $ print_help()
  stop()
} else {
  pargs <- parser_qc $ parse_args(vargs)
}

print(parser_qc $ command)
