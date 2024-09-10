
parser <- argparse::ArgumentParser(
  prog = "run",
  description = "run custom script"
)

parser $ add_argument(
  type = "character", dest = "file",
  help = "path to the script file"
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

if (!file.exists(pargs $ file)) {
  cat(red(paste("file", pargs $ file, "does not exist")), crlf)
  stop()
}

source(pargs $ file)
