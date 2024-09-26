
parser <- argparse::ArgumentParser(
  prog = "source",
  description = "run commands recorded in script"
)

parser $ add_argument(
  "-s", type = "character", dest = "script", default = "",
  help = "the input script file, will be executed line by line"
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

cmdlist <- c(cmdlist, parse_script(pargs $ script))
