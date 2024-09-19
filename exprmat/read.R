
if (shared[["is_loaded"]]) {
  cat("data are already loaded in this directory.", crlf)
  stop()
}

if (!shared[["is_reference_assigned"]]) {
  cat("you should assign a reference genome first.", crlf)
  stop()
}

cat(stringr::str_wrap(indent = 0, exdent = 0, width = 80, paste(
  "this module processed an regular input of a single merged rna-seq data from",
  "a single expression matrix and necessary sample descriptions and gene",
  "descriptions. all these information should be prepared to a tab delimited",
  "text file in advance. the expression matrix only receives raw counts. the",
  "normalization and quality control process will be handled here. for",
  "multiple datasets, you should merge them (e.g. prepare qc for each of them,",
  "finding variable genes, finding anchors, integrating anchors, and extract",
  "only the merged dataset to a raw-count expression matrix.)"
)), crlf, crlf)

cat(stringr::str_wrap(indent = 0, exdent = 0, width = 80, paste(
  magenta("expression matrix in raw counts:"), "specify an expression matrix",
  "file, text files or compressed text files are accepted. should contain a",
  yellow("column header"), "with tabs, and comment lines start with '#'"
)), crlf)

fname_expr_count <- read()

expr_count <- read.delim(fname_expr_count, sep = "\t", header = TRUE,
                         comment.char = "#")

expr_count_rownames <- rownames(expr_count)
expr_count <- tibble(expr_count)

cat(crlf)
if (!is.null(expr_count_rownames)) {
  cat("row names detected: of size", length(expr_count_rownames), crlf)
  print(head(expr_count_rownames))
}

if (!is.null(expr_count |> colnames())) {
  cat("column names:", length(expr_count_rownames), crlf)
  print(head(expr_count |> colnames()))
}

cat(crlf)
print(expr_count)

cat(crlf)

cat(stringr::str_wrap(indent = 0, exdent = 0, width = 80, paste(
  magenta("gene mapping:"), "the following questions will request your",
  "configuration on gene mapping, and we will automatically assign your genes",
  "to the predefined database by refseq and ensembl."
)), crlf, crlf)

cat(stringr::str_wrap(indent = 0, exdent = 4, width = 80, paste(
  magenta("1) "), "which columns are your gene metadata?"
)), crlf, crlf)

gene_meta_cols <- read()
gene_meta_cols <- str_split(gene_meta_cols, " ")[[1]]

for (cx in gene_meta_cols) {
  if (!(cx %in% colnames(expr_count))) {
    cat(crlf, red("error:"), "the specified column", c, "does not exist!", crlf)
    q(save = "no", status = 1)
  }
}

cat(crlf)

cat(
  stringr::str_wrap(indent = 0, exdent = 4, width = 80, paste(
    magenta("2) "), "which type of pivot column you would like to use when",
    "assigning your gene to the database? possible selections:"
  )),
  "\n    (a) 'entrez': ncbi entrez gene id, integers",
  "\n    (b) 'hgnc': the refseq name, full capitalized",
  "\n    (c) 'ensembl': the ensembl gene id, starting with ENSG or ENSMUSG",
  crlf, crlf
)

gene_map_method <- read()

cat(crlf)

cat(stringr::str_wrap(indent = 0, exdent = 4, width = 80, paste(
  magenta("3) "), "which column is the pivot column?"
)), crlf, crlf)

gene_map_pivot <- read()

if (!(gene_map_pivot %in% gene_meta_cols)) {
  cat(crlf, red("error:"), "the specified pivot column", gene_map_pivot,
      "doesn't exist!", crlf)
  q(save = "no", status = 1)
}

pivot_name <- gene_map_pivot
gene_map_pivot <- pull(expr_count, gene_map_pivot)

if (is.na(gene_map_pivot) |> sum() > 0)
  cat(crlf, yellow("the pivot column contains NAs!"), crlf,
      "   contains", red(is.na(gene_map_pivot) |> sum()), "NAs", crlf)

genes <- readRDS("genome.rds")
gene_names <- genes $ name

genes_notfound <- c()
genes_duplicate <- c()
genes_mask <- c()
genes_map <- c()

if (gene_map_method == "entrez") {

  # the entrez ids are the original primary key in out gene annotation table.
  # the entrez id are not duplicated. so it is the most striaght-forward way.

  m <- match(gene_map_pivot, genes $ entrez)

  # if the gene mapping is not found in reference, a NA will occur.
  # if duplicate match, the same index is called again.
  # if the reference is duplicated, will only return the first match.

  genes_notfound <- gene_map_pivot[is.na(m)]
  genes_duplicate <- gene_map_pivot[duplicated(m[!is.na(m)])]
  genes_mask <- !is.na(m)
  genes_map <- m[!is.na(m)]

  expr_count <- expr_count[genes_mask, ]
  genes_meta <- genes[genes_map, ]

  expr_count <- expr_count[!duplicated(genes_meta $ entrez), ]
  genes_meta <- genes_meta[!duplicated(genes_meta $ entrez), ]

  # may have duplicated names (tRNAs, see below), but no duplicated entrez.
  # for these genes, we will just keep the first one.

  # HOWEVER, YOU SHOULD MANUALLY CHECK WHETHER THE DUPLICATED NAMES ARE ONE
  # THING. this get better resolved in the recent releases of annotations,
  # however, some annotations may use the same name as completely different
  # abbreviations of gene full name.

  # "
  #     However, I can say definitively that you should NOT merge the
  #     expression values of genes that have the same symbol. Gene symbols
  #     seem to shuffle around more haphazardly than Ensembl IDs or Entrez
  #     IDs, and so it is pretty common for completely different genes to
  #     be labeled with the same gene symbol depending on which gene
  #     annotation builds are being used.
  #
  #     Which brings me to my next point. You will save yourself some
  #     headaches if you know what gene annotation build was used in
  #     originally mapping your RNA-seq reads to Ensembl IDs. Use the
  #     same build to translate from Ensembl IDs to gene symbols. This
  #     won't necessarily get rid of all issues, but it will be a big help.
  #                                                                          "
  #                      <https://www.researchgate.net/post/How-to-deal-with-
  #                       multiple-ensemble-IDs-mapping-to-one-gene-symbol-in-
  #                       a-RNA-Seq-dataset>

} else if (gene_map_method == "hgnc") {

  # in this case, we will face disturbing naming problems. some of the genes
  # do not have the same name. it should be noted that mitochondrial genes are
  # especially tend to have others names.

  # in the refseq database, mitochondrial gene names are

  mitos_ncbi <- c(
    "COX1", "COX2", "COX3", "ND1", "ND2", "ND3", "ND4L", "ND4", "ND5", "ND6",
    "CYTB", "ATP6", "ATP8"
  )

  mitos_2 <- paste("MT-", mitos_ncbi, sep = "")

  # don't be confused with these cox1/2/3. they represent cytochrome oxidase.
  # the inflammatory cytokine synthase cox-1/2 are officially named
  # prostoglandin endoperoxide synthase 1/2 or ptgs1/2 actually ...

  # for some occasion, they are

  mitos_3 <- c(
    "mt-Co1", "mt-Co2", "mt-Co3",
    "mt-Nd1", "mt-Nd2", "mt-Nd3", "mt-Nd4l", "mt-Nd4", "mt-Nd5", "mt-Nd6",
    "mt-Cytb", "mt-Atp6", "mt-Atp8"
  )

  # where cytochrome-c oxidase can be named mt-co1/2/3 but not cox1/2/3. anyway,
  # these issues suggest that we should use indices as much as possible, since
  # the names are constantly varying. you may encounter false misses when using
  # names to specify genes.

  # we do not use the former lists as mitochondrial gene names. these are only
  # shown as examples of what contained in the mt genome. in our new reference
  # construction procedures, we has a `mito` column indicating whether the gene
  # is a mitochondrial gene. (do not consider plants though ;-) ...

  # for human there are roughly five cases (or two genes, both tRNA-encoding)
  # faces the condition that they have shared name but with different entrez.
  # the mouse genome have not even one. i think we can disgard them.

  genes_chr <- genes[!genes $ mito, ]

  m <- match(gene_map_pivot, genes_chr $ gene)
  genes_notfound <- gene_map_pivot[is.na(m)] # mito genes should not be found.
  genes_duplicate <- gene_map_pivot[duplicated(m[!is.na(m)])]
  genes_mask <- !is.na(m)
  genes_map <- m[!is.na(m)]

  expr_count_chr <- expr_count[genes_mask, ]
  genes_meta <- genes_chr[genes_map, ]

  expr_count_chr <- expr_count_chr[!duplicated(genes_meta $ entrez), ]
  genes_meta <- genes_meta[!duplicated(genes_meta $ entrez), ]

  # process mitochondrial genes separately

  genes_mito <- genes[genes $ mito, ]
  mmt <- match(gene_map_pivot, genes_mito $ gene)
  if (sum(!is.na(mmt)) == 0) # all the genes not found. suggesting a name error
    mmt <- match(gene_map_pivot, paste("mt-", genes_mito $ gene))
  if (sum(!is.na(mmt)) == 0) # again
    mmt <- match(gene_map_pivot, paste("MT-", genes_mito $ gene))
  if (sum(!is.na(mmt)) == 0) # again
    mmt <- match(gene_map_pivot, paste("Mt-", genes_mito $ gene))
  if (sum(!is.na(mmt)) == 0) # again
    mmt <- match(gene_map_pivot, paste("MT", genes_mito $ gene))
  if (sum(!is.na(mmt)) == 0) # again
    mmt <- match(gene_map_pivot, paste("mt", genes_mito $ gene))
  if (sum(!is.na(mmt)) == 0) # again
    mmt <- match(gene_map_pivot, paste("Mt", genes_mito $ gene))

  mt_mask <- !is.na(mmt)
  mt_map <- mmt[!is.na(mmt)]

  expr_count_mt <- expr_count[mt_mask, ]
  genes_meta_mt <- genes_mito[mt_map, ]

  expr_count_mt <- expr_count_mt[!duplicated(genes_meta_mt $ entrez), ]
  genes_meta_mt <- genes_meta_mt[!duplicated(genes_meta_mt $ entrez), ]

  newfound_name <- expr_count_mt |> pull(pivot_name)

  cat(crlf, yellow("found these mitochondrial genes:"), crlf)
  if (length(newfound_name) > 0) {
    cat("    ")
    print(newfound_name)
  } else {
    cat("   ", red("not a single mitochondrial genes found, check your data!"))
    cat(crlf)
  }

  genes_notfound <- base::setdiff(genes_notfound, newfound_name)
  expr_count <- rbind(expr_count_mt, expr_count_chr)
  genes_meta <- rbind(genes_meta_mt, genes_meta)

} else if (gene_map_method == "ensembl") {

  dcensem <- decollapse(genes, "ensembl", sep = ";")
  m <- match(gene_map_pivot, dcensem $ ensembl)

  # if the gene mapping is not found in reference, a NA will occur.
  # if duplicate match, the same index is called again.
  # if the reference is duplicated, will only return the first match.

  genes_notfound <- gene_map_pivot[is.na(m)]
  genes_duplicate <- gene_map_pivot[duplicated(m[!is.na(m)])]
  genes_mask <- !is.na(m)
  genes_map <- m[!is.na(m)]

  expr_count <- expr_count[genes_mask, ]
  genes_meta <- dcensem[genes_map, ]

  expr_count <- expr_count[!duplicated(genes_meta $ ensembl), ]
  genes_meta <- genes_meta[!duplicated(genes_meta $ ensembl), ]

} else {
  cat(crlf, red("error:"), "invalid mapping method!", crlf)
  q(save = "no", status = 1)
}

# print duplicated names
dup_name <- genes_meta[duplicated(genes_meta $ gene), ] $ gene |> unique()
dup_name <- genes_meta[
  pull(genes_meta, "gene") %in% dup_name,
  c("gene", "entrez", ".start", ".end", ".strand", ".seqid")
]

if (nrow(dup_name) > 0) {
  cat(crlf, red("these gene names got duplicated!"), crlf, crlf)
  print(dup_name)
}

# print mapping result.
cat(crlf, expr_count |> nrow() |> green() |> italic(),
    italic("genes assigned."), crlf)

if (genes_notfound |> length() > 0)
  cat(
    crlf, genes_notfound |> length() |> red() |> italic(),
    italic("genes detected but not found in reference genome.\n   "),
    genes_notfound |> head() |> str_c(collapse = " ") |> red() |> italic(),
    "...", crlf
  )

if (genes_duplicate |> length() > 0)
  cat(
    crlf, genes_duplicate |> length() |> red() |> italic(),
    italic("genes duplicated.\n   "),
    genes_duplicate |> head() |> str_c(collapse = " ") |> red() |> italic(),
    "...", crlf
  )

for (cx in gene_meta_cols) {
  genes_meta[, paste("user_", cx, sep = "")] <- expr_count[, cx]
}

if (nrow(genes_meta) != nrow(expr_count)) {
  cat(crlf, red("error:"), "assertion failed, genes_meta has", nrow(genes_meta),
      "columns, while expr_count has", nrow(expr_count), "columns", crlf)
  q(save = "no", status = 1)
}

fname_sample_meta <- "sample.tsv"

# read the sample data. the id row is required, and must match the header.
sample_meta <- read.delim(fname_sample_meta, sep = "\t", header = TRUE,
                          comment.char = "#")

n_sample <- nrow(sample_meta)
columns <- colnames(expr_count)
sample_columns <- columns[!(columns %in% gene_meta_cols)]
sample_columns <- sample_columns[!duplicated(sample_columns)]

expr_count <- expr_count[, sample_columns]

if (length(sample_columns) != n_sample) {
  cat(crlf, red("error:"), "mismatch between sample config",
      "and real sample columns.", length(sample_columns), "!=", n_sample, crlf)
  q(save = "no", status = 1)
}

if (setdiff(sample_columns, sample_meta $ id) |> length() > 0) {
  cat(crlf, red("error:"), "sample indexes mismatch", crlf)
  q(save = "no", status = 1)
}

ord <- c()
for (cx in sample_columns) {
  id <- which(sample_meta $ id == cx)
  ord <- c(ord, id[1])
}

sample_meta <- sample_meta[ord, ]

# display the final processed data that is capable to form an adequate
# single cell experiment expression matrix.

dir.create("features")

sample_meta <- sample_meta |> tibble()
genes_meta <- genes_meta |> tibble()

saveRDS(sample_meta, "features/samples-meta.rds")
saveRDS(genes_meta, "features/genes-meta.rds")
saveRDS(expr_count, "features/matrix.rds")

cat(crlf, green("successfully saved expression matrix."), crlf)

# remove private variables generated in the script.
suppressWarnings(rm(
  fname_expr_count, expr_count_rownames, gene_meta_cols, gene_map_method,
  gene_map_pivot, gene_names, genes_duplicate, genes_map, genes_mask,
  genes_missing, genes_notfound, ids, mitos_2, mitos_3, mitos_id, mitos_ncbi,
  fname_sample_meta, n_sample, columns, sample_columns, id, ord, genes
))

shared[["counts"]] <- expr_count
shared[["meta_sample_raw"]] <- sample_meta
shared[["meta_gene_raw"]] <- genes_meta
