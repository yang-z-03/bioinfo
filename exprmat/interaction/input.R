
parser <- argparse::ArgumentParser(
  prog = "cpdbin",
  description = "generate the inputs required by cellphonedb"
)

parser $ add_argument(
  "--subset", dest = "cond",
  type = "character", default = c(), nargs = "+",
  help = paste("conditions specifying the cluster and groupings")
)

parser $ add_argument(
  "-c", dest = "cluster", type = "character", default = "",
  help = paste("specify the cluster identification")
)

parser $ add_argument(
  "-d", dest = "de", type = "character", default = "",
  help = paste("specify the differentially expressed gene names (a key from",
               "the shared object table)")
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

auto_cast <- function(binary) {
  destclass <- shared[["seurat"]] @ meta.data |>
    dplyr::pull(binary[1]) |>
    class()

  if (destclass == "factor")
    destclass <- shared[["seurat"]] @ meta.data |>
      dplyr::pull(binary[1]) |>
      levels() |>
      class()

  if (destclass == "character") {
    comp <- as.character(binary[2])
  } else if (destclass == "logical") {
    comp <- as.logical(binary[2])
  } else {
    comp <- as.numeric(binary[2])
  }

  return(comp)
}

parse_subset_cond <- function(cond) {
  if (stringr::str_detect(cond, ">=")) {
    binary <- stringr::str_split_1(cond, ">=")
    comp <- auto_cast(binary)
    return(shared[["seurat"]][[binary[1]]] >= comp)

  } else if (stringr::str_detect(cond, "<=")) {
    binary <- stringr::str_split_1(cond, "<=")
    comp <- auto_cast(binary)
    return(shared[["seurat"]][[binary[1]]] <= comp)

  } else if (stringr::str_detect(cond, ">")) {
    binary <- stringr::str_split_1(cond, ">")
    comp <- auto_cast(binary)
    return(shared[["seurat"]][[binary[1]]] > comp)

  } else if (stringr::str_detect(cond, "<")) {
    binary <- stringr::str_split_1(cond, "<")
    comp <- auto_cast(binary)
    return(shared[["seurat"]][[binary[1]]] < comp)

  } else if (stringr::str_detect(cond, "=")) {
    binary <- stringr::str_split_1(cond, "=")
    comp <- auto_cast(binary)
    return(shared[["seurat"]][[binary[1]]] == comp)
  }
}

parse_subset_single <- function(arg) {
  is_combined <- stringr::str_detect(arg, "&")
  if (is_combined) {
    splits <- stringr::str_split_1(arg, "&")
    condlist <- lapply(splits, parse_subset_cond)
    return(purrr::reduce(condlist, `&`))
  } else {
    return(parse_subset_cond(arg))
  }
}

parse_subset <- function(arg) {
  condlist <- lapply(arg, parse_subset_single)
  return(purrr::reduce(condlist, `|`))
}

if (pargs $ cond |> length() == 0) {
  filter <- rep(TRUE, ncol(shared[["seurat"]]))
} else {
  filter <- parse_subset(pargs $ cond)
}

filter_cell <- shared[["meta_sample"]][filter, ]
cellnames <- SeuratObject::Cells(filter_expr)

shared[["seurat"]] $ .temp.filter <- filter
filter_expr <- shared[["seurat"]] |> subset(subset = .temp.filter == TRUE)
shared[["seurat"]] $ .temp.filter <- NULL

# seurat requires a non-unique gene name for each feature. while this should
# also be compatible with hgnc names, the few (3 in human genome) genes whose
# name contains a underscore (automatically turned away by seurat) will be lost.
# however, i think i do not need those minorities :)
# the cell names are auto-generated to be c1 c2 ... etc.

#. counts <- SeuratObject::GetAssayData(filter_expr, layer = "counts")

# cellphone db requires normalized data in linear space.
# we can not determine whether the data stored in 'counts' are normalized or not
# so we just exponentiate the data layer.

norm <- SeuratObject::GetAssayData(filter_expr, layer = "data") |> exp() |> t()
norm <- norm |> as("sparseMatrix")
rownames(norm) <- paste("c", seq_along(cellnames), sep = "")
colnames(norm) <- shared[["meta_gene"]] $ seurat_names

# for anndata, observations (rows) are cells, and variables (columns) are genes
# so we should transpose the seurat matrix.

require(anndata)

gene_ensembl_expand <- decollapse(shared[["meta_gene"]], "ensembl", ";")
gene_select <- match(
  shared[["meta_gene"]] $ gene,
  gene_ensembl_expand $ gene
)
ensembl_list <- gene_ensembl_expand[gene_select, ] $ ensembl

ann <- anndata::AnnData(
  X = norm |> as("sparseMatrix"),
  obs = data.frame(
    name = cellnames,
    cluster = pull(filter_cell, pargs $ cluster),
    row.names = paste("c", seq_along(cellnames), sep = "")
  ),
  var = data.frame(
    seurat = shared[["meta_gene"]] $ seurat_names,
    hgnc = shared[["meta_gene"]] $ gene,
    entrez = shared[["meta_gene"]] $ entrez,
    ensembl = ensembl_list,
    row.names = shared[["meta_gene"]] $ seurat_names
  )
)

meta <- data.frame(
  barcode_sample = paste("c", seq_along(cellnames), sep = ""),
  cell_type = pull(filter_cell, pargs $ cluster)
)

keys <- names(shared[[pargs $ de]])
ctype <- c()
de <- c()
for (kname in keys) {
  ctype <- c(ctype, rep(kname, shared[[pargs $ de]][[kname]] |> length()))
  de <- c(de, shared[[pargs $ de]][[kname]])
}

deg <- data.frame(
  cell_type = ctype,
  gene = de
)

# save the files

if (!dir.exists("cpdb")) dir.create("cpdb")

cat(blue("writing metadata table ..."), crlf)
write.table(meta, file = "cpdb/metadata.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

cat(blue("writing differential expressing genes table ..."), crlf)
write.table(deg, file = "cpdb/deg.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

cat(blue("writing annotated data ..."), crlf)
suppressWarnings(suppressMessages(
  write_h5ad(ann, "cpdb/norm-counts.h5ad")
))
