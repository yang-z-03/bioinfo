
parser <- argparse::ArgumentParser(
  prog = "annot",
  description = "annotate cell types"
)

parser $ add_argument(
  "-c", dest = "cluster", type = "character", default = "",
  help = paste("specify the cluster identification. you may also skip this",
               "option to classify every cell separately")
)

parser $ add_argument(
  "-r", dest = "refs",
  type = "character", default = c(), nargs = "+",
  help = paste("names to the reference datasets (harmonized, cell-ontology",
               "labeled, single cell experiment or bulk experiment). thers",
               "datasets are stored in ctref/* directory, where they are",
               "generated and format using 'ctref'")
)

parser $ add_argument(
  "--pseudo-bulk", dest = "psbulk", action = "store_true", default = FALSE,
  help = paste("make pseudo-bulk references")
)
