
parser <- argparse::ArgumentParser(
  prog = "gffselect",
  description = "select subset columns and types from original gff"
)

parser $ add_argument(
  type = "character", dest = "taxo",
  help = "generate gff subset files for given taxo"
)

parser $ add_argument(
  "-t", dest = "type", type = "character", nargs = "*", default = c(),
  help = "selected annotation types"
)

parser $ add_argument(
  "-c", dest = "columns", type = "character", nargs = "*", default = c(),
  help = "selected columns from this type, by default dumps all"
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

cmdlineargs <- paste(
  paste(gp_annot, "gffselect", sep = "/"),
  "-z", paste(gp_refseq, pargs $ taxo, "genomic.gff.gz", sep = "/"),
  "-t", stringr::str_c(pargs $ type, collapse = " "),
  "-c", stringr::str_c(pargs $ columns, collapse = " ")
)

cat(cmdlineargs, crlf)

system(cmdlineargs)
