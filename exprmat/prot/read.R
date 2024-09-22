
parser <- argparse::ArgumentParser(
  prog = "readp",
  description = "read quantitative proteomic data"
)

parser $ add_argument(
  "-n", dest = "name", type = "character", default = "uniprot",
  help = "key for protein identifer. expect 'uniprot'"
)

parser $ add_argument(
  "-f", dest = "file", type = "character", default = "filename",
  help = "the quantitative matrix, in tab-delimited format"
)

parser $ add_argument(
  "-t", dest = "trans",
  default = FALSE, action = "store_true",
  help = paste("transpose the matrix (the input matrix should have proteins as",
               "rows, and samples as columns, if not, you should transpose it")
)

parser $ add_argument(
  "--prot-meta", dest = "meta.p", type = "character",
  default = c(), nargs = "*",
  help = "extract protein (row) metadata from columns"
)

parser $ add_argument(
  "--sample-meta", dest = "meta.s", type = "character",
  default = c(), nargs = "*",
  help = "extract sample (column) metadata from rows"
)

parser $ add_argument(
  "--has-pname", dest = "has.p",
  default = FALSE, action = "store_true",
  help = paste("the table has protein names (as the first column)")
)

parser $ add_argument(
  "--has-sname", dest = "has.s",
  default = FALSE, action = "store_true",
  help = paste("the table has sample names (as the first row)")
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

# the purpose is to read an un-normalized protein quantification data
# and store in the similar seurat/protein-info/sample-info mode.
# that is ready for later gsea or de analyses.

table <- read.table(
  pargs $ file, header = FALSE, sep = "\t",
  numerals = "allow.loss", encoding = "UTF-8"
)

if (pargs $ trans) {
  table <- table |> data.table::transpose()
}

if (pargs $ has.p) {
  rownames(table) <- table $ V1
} else {
  rownames(table) <- paste("r", seq_len(nrow(table)), sep = "")
}

if (pargs $ has.s) {
  colnames(table) <- table[1, ]
} else {
  colnames(table) <- paste("c", seq_len(ncol(table)), sep = "")
}

# now extract the metadata rows and columns

meta_sample <- data.frame(
  pivot = colnames(table),
  .id = colnames(table)
)

if (pargs $ meta.s |> length() > 0) {
  for (x in pargs $ meta.s) {
    meta_sample[[x]] <- table[x, ] |> unlist()
  }
}

meta_protein <- data.frame(
  pivot = rownames(table),
  .id = rownames(table)
)

if (pargs $ meta.p |> length() > 0) {
  for (x in pargs $ meta.p) {
    meta_protein[[x]] <- table |> dplyr::pull(x)
  }
}

remove_rows <- length(pargs $ meta.s) + if (pargs $ has.s) 1 else 0
remove_cols <- length(pargs $ meta.p) + if (pargs $ has.p) 1 else 0

meta_protein <- meta_protein[-(1:remove_rows), ] |> tibble()
meta_sample <- meta_sample[-(1:remove_cols), ] |> tibble()
table <- table[-(1:remove_rows), -(1:remove_cols)] |> tibble()

# mapping to the genome annotation

gt <- readRDS("genome.rds")

switch(
  pargs $ name,
  uniprot = {

    cat(blue("duplicates in the querying uniprot id:"), crlf)
    print(meta_protein $ .id |> duplicated() |> table())
    cat(crlf)
    dupmask <- meta_protein $ .id |> duplicated()
    meta_protein <- meta_protein[!dupmask, ]
    table <- table[!dupmask, ]

    gtprot <- gt[, c("entrez", "uniprot")]
    gtprot <- decollapse(gtprot, "uniprot", ";")

    matched <- match(meta_protein $ .id, gtprot $ uniprot)
    cat(blue("not found in the genomic annotation:"), crlf)
    print(matched |> is.na() |> table())

    entrez_id <- gtprot[matched, ] $ entrez
    gt_result <- match(entrez_id, gt $ entrez)
    gt_result <- gt[gt_result, ]

    meta_protein <- cbind(meta_protein, gt_result) |> tibble()
    naprots <- is.na(meta_protein $ entrez)
    meta_protein <- meta_protein[!naprots, ]
    table <- table[!naprots, ]
  }
)

shared[["proteins"]] <- table
shared[["meta_sample"]] <- meta_sample
shared[["meta_protein"]] <- meta_protein

if (!dir.exists("proteome")) dir.create("proteome")

saveRDS(table, "proteome/counts.rds")
saveRDS(meta_sample, "proteome/samples-meta.rds")
saveRDS(meta_protein, "proteome/proteins-meta.rds")
