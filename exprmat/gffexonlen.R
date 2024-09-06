
parser <- argparse::ArgumentParser(
  prog = "gffexonlen",
  description = "count exon lengths"
)

parser $ add_argument(
  type = "character", dest = "taxo",
  help = "generate gff index files for given taxo"
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

cmdlineargs <- paste(
  paste(gp_annot, "gffexonlen", sep = "/"),
  "-z", paste(gp_refseq, pargs $ taxo, "genomic.gff.gz", sep = "/")
)

system(cmdlineargs)
