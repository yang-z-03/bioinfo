
parser_sc3 <- argparse::ArgumentParser(prog = "cluster")

parser_sc3 $ add_argument(type = "character", dest = "command",
                          help = "clustering commands")

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser_sc3 $ print_help()
  stop()
} else {
  pargs <- parser_sc3 $ parse_args(vargs)
}

print(pargs $ command)
