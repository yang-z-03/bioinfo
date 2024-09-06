
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

cat(stringr::str_wrap(indent = 0, exdent = 4, width = 80, paste(
  magenta("2) "), "which type of pivot column you would like to use when",
  "assigning your gene to the database? possible selections:",
  "\n(a) 'entrez': ncbi entrez gene id, integers;",
  "\n(b) 'refseq' or 'capitalize': the refseq sequence name, full capitalized;",
  "\n(c) 'ensembl': the ensembl gene id, starting with ENSG or ENSMUSG;"
)), crlf, crlf)

#  [1] "entrez"             "refseq_seqname"     "refseq_source"
#  [4] "refseq_feature"     "refseq_start"       "refseq_end"
#  [7] "refseq_score"       "refseq_strand"      "refseq_frame"
# [10] "refseq_id"          "refseq_xref"        "refseq_name"
# [13] "refseq_description" "refseq_genbank"     "refseq_gene"
# [16] "refseq_biotype"     "genbank_accession"  "symbol_alias"
# [19] "ensembl"            "ensembl_prot"       "ensembl_transcript"
# [22] "ec"                 "go_evidence"        "go_evidence_all"
# [25] "name"               "type"               "go"
# [28] "go_all"             "ipi_accession"      "map_location"
# [31] "omim"               "go_ontology"        "go_ontology_all"
# [34] "kegg_pathway"       "pfam"               "pmid"
# [37] "prosite"            "refseq"             "symbol"
# [40] "ucsc_known_gene"    "uniprot"

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

gene_map_pivot <- pull(expr_count, gene_map_pivot)
gene_names <- genes $ refseq_name
genes_notfound <- c()
genes_missing <- c()
genes_map <- c() # mapping order to genes
genes_duplicate <- c()
genes_mask <- c() # mapping order to pivot

if (gene_map_method == "entrez") {

  for (t in gene_map_pivot) {

    if (is.na(t)) {
      genes_mask <- c(genes_mask, FALSE)
      next
    }

    ids <- grep(paste("^", t, "$", sep = ""),
                genes $ entrez, ignore.case = FALSE)
    if (ids |> length() == 0) {
      genes_notfound <- c(genes_notfound, t)
      genes_mask <- c(genes_mask, FALSE)
    } else if (ids |> length() == 1) {
      genes_map <- c(genes_map, ids)
      genes_mask <- c(genes_mask, TRUE)
    } else {
      genes_map <- c(genes_map, ids[1])
      genes_duplicate <- c(genes_duplicate, ids[2:length(ids)])
      genes_mask <- c(genes_mask, TRUE)
    }
  }
  genes_missing <- setdiff(seq_along(genes), c(genes_map, genes_duplicate))

} else if (gene_map_method == "refseq" || gene_map_method == "capitalize") {

  # in this case, we will face disturbing naming problems. some of the genes
  # do not have the same name. it should be noted that mitochondrial genes are
  # especially tend to have others names.

  # in the refseq database, mitochondrial gene names are

  mitos_ncbi <- c(
    "COX1", "COX2", "COX3", "ND1", "ND2", "ND3", "ND4L", "ND4", "ND5", "ND6",
    "CYTB", "ATP6", "ATP8"
  )

  mitos_2 <- paste("mt-", mitos_ncbi, sep = "")

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

  mitos_id <- c()
  for (m in mitos_ncbi) {
    mitos_id <- c(mitos_id, grep(paste("^", m, "$", sep = ""),
                                 genes $ refseq_name))
  }

  for (t in gene_map_pivot) {

    if (is.na(t)) {
      genes_mask <- c(genes_mask, FALSE)
      next
    }

    ids <- grep(paste("^", t, "$", sep = ""),
                genes $ refseq_name, ignore.case = TRUE)

    if (startsWith(str_to_lower(t), "mt")) {
      ids <- mitos_id[
        grep(paste("^", t, "$", sep = ""), mitos_2, ignore.case = TRUE)
      ]

      if (length(ids) == 0) ids <- mitos_id[
        grep(paste("^", t, "$", sep = ""), mitos_3, ignore.case = TRUE)
      ]
    }

    if (ids |> length() == 0) {
      genes_notfound <- c(genes_notfound, t)
      genes_mask <- c(genes_mask, FALSE)
    } else if (ids |> length() == 1) {
      genes_map <- c(genes_map, ids)
      genes_mask <- c(genes_mask, TRUE)
    } else {
      genes_map <- c(genes_map, ids[1])
      genes_duplicate <- c(genes_duplicate, ids[2:length(ids)])
      genes_mask <- c(genes_mask, TRUE)
    }
  }
  genes_missing <- setdiff(seq_along(genes), c(genes_map, genes_duplicate))

} else if (gene_map_method == "ensembl") {

  for (t in gene_map_pivot) {

    if (is.na(t)) {
      genes_mask <- c(genes_mask, FALSE)
      next
    }

    ids <- grep(paste("^", t, "$", sep = ""),
                genes $ ensembl, ignore.case = TRUE)
    if (ids |> length() == 0) {
      genes_notfound <- c(genes_notfound, t)
      genes_mask <- c(genes_mask, FALSE)
    } else if (ids |> length() == 1) {
      genes_map <- c(genes_map, ids)
      genes_mask <- c(genes_mask, TRUE)
    } else {
      genes_map <- c(genes_map, ids[1])
      genes_duplicate <- c(genes_duplicate, ids[2:length(ids)])
      genes_mask <- c(genes_mask, TRUE)
    }
  }
  genes_missing <- setdiff(seq_along(genes), c(genes_map, genes_duplicate))

} else {
  cat(crlf, red("error:"), "invalid mapping method!", crlf)
  q(save = "no", status = 1)
}

# print mapping result.
cat(crlf, genes_map |> length() |> green() |> italic(),
    italic("genes assigned."), crlf)

if (genes_notfound |> length() > 0)
  cat(
    crlf, genes_notfound |> length() |> red() |> italic(),
    italic("genes detected but not found in reference genome.\n   "),
    genes_notfound |> head() |> str_c(collapse = " ") |> red() |> italic(),
    "...", crlf
  )

if (genes_missing |> length() > 0)
  cat(
    crlf, genes_missing |> length() |> yellow() |> italic(),
    italic("genes not detected.\n   "),
    gene_names[genes_missing |> head()] |> str_c(collapse = " ") |> yellow() |> italic(), # nolint
    "...", crlf
  )

if (genes_duplicate |> length() > 0)
  cat(
    crlf, genes_duplicate |> length() |> cyan() |> italic(),
    italic("duplicated genes in the database.\n   "),
    gene_names[genes_duplicate |> head()] |> str_c(collapse = " ") |> cyan() |> italic(), # nolint
    "...", crlf
  )

expr_count <- expr_count[genes_mask, ]
genes_meta <- genes[genes_map, ]

for (cx in gene_meta_cols) {
  genes_meta[, paste("user_", cx, sep = "")] <- expr_count[, cx]
}

if (nrow(genes_meta) != nrow(expr_count)) {
  cat(crlf, red("error:"), "assertion failed, genes_meta has", nrow(genes_meta),
      "columns, while expr_count has", nrow(expr_count), "columns", crlf)
  q(save = "no", status = 1)
}

cat(crlf)

cat(stringr::str_wrap(indent = 0, exdent = 0, width = 80, paste(
  magenta("sample table:"), "the sample and experiment configuration.",
  "will be appended to the column metadata. a tab-delimited text is expected,",
  "and a column named", yellow("'id'"), "is required and must match the header",
  "of the expression matrix."
)), crlf)

fname_sample_meta <- read()

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

cat(crlf, green("successfully saved expression matrix."), crlf, crlf)

# remove private variables generated in the script.
suppressWarnings(rm(
  fname_expr_count, expr_count_rownames, gene_meta_cols, gene_map_method,
  gene_map_pivot, gene_names, genes_duplicate, genes_map, genes_mask,
  genes_missing, genes_notfound, ids, mitos_2, mitos_3, mitos_id, mitos_ncbi,
  fname_sample_meta, n_sample, columns, sample_columns, id, ord
))

gc()
