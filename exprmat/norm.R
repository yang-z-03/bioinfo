
parser <- argparse::ArgumentParser(
  prog = "norm",
  description = "normalization and transformation of counts"
)

parser $ add_argument(
  "-m", dest = "method", type = "character",
  help = paste("method of normalization. supports 'cpm', 'fpkm', 'tpm', 'uq',",
               "'rle', 'tmm', 'seurat', and 'sct'")
)

parser $ add_argument(
  "-r", dest = "regress", type = "character", nargs = "*", default = NULL,
  help = paste("select which columns of confounders to regress out. primarily",
               "quality control statistics such as 'reads' (the total read",
               "depth for each sample), 'detection', 'pct_ribo', 'pct_mito',",
               "'doublet_scores', 's_score' and 'g2m_score'",
               "(estimated cell cycle status)")
)

parser $ add_argument(
  "--p", dest = "uq", type = "double", default = 0.99,
  help = paste("p-value cutoff for upper quartile. (by default 0.99)")
)

parser $ add_argument(
  "--trim-m", dest = "m", type = "double", default = 0.3,
  help = paste("trim to use on log-ratios for tmm.")
)

parser $ add_argument(
  "--trim-a", dest = "a", type = "double", default = 0.05,
  help = paste("trim to use on the combined absolute levels for tmm.")
)

parser $ add_argument(
  "--cutoff-a", dest = "ca", type = "double", default = -1e10,
  help = paste("cutoff on A values before trimming for tmm")
)

parser $ add_argument(
  "--no-weight", dest = "weighting", action = "store_false", default = TRUE,
  help = paste("compute (asymptotic binomial precision) weights for tmm")
)

parser $ add_argument(
  "--no-cc", dest = "cc", action = "store_false", default = TRUE,
  help = paste("compute cell cycle scores and enable such regresses")
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

if (!shared[["is_qc"]]) {
  cat("you should run 'qc' first", crlf)
  stop()
}

taxo <- readRDS("taxo.rds")

if (length(taxo $ g2m) == 0 || length(taxo $ s) == 0) {
  if (pargs $ cc) {
    cat(str_wrap(
      yellow("your reference genome doesn't support cell cycle scoring."),
      "you should rerun refer command with --s-phase and --g2m-phase.",
      "we have turned on --no-cc automatically for you."
    ))
  }
  pargs $ cc <- FALSE
}

if (!pargs $ cc && !is.null(pargs $ regress)) {
  if ("s_score" %in% pargs $ regress || "g2m_score" %in% pargs $ regress) {
    cat("cell cycle scoring is not available for current reference.", crlf)
    cat("you should rerun refer with --s-phase and --g2m-phase.", crlf)
    stop()
  }
}

expr_count <- readRDS("qc/matrix.rds")
expr_mat <- expr_count |> as("sparseMatrix")
sample_meta <- readRDS("qc/samples-meta.rds")
genes_meta <- readRDS("qc/genes-meta.rds")
exonlen <- read.delim("exonlens.tsv")

genes_meta $ .key <-
  paste(genes_meta $ ensembl, genes_meta $ ensembl_version, sep = ".")

genes_meta <- merge(
  x = genes_meta, y = exonlen[, c("symbol", "mean", "median", "merged")],
  by.x = ".key", by.y = "symbol",
  all.x = TRUE, all.y = FALSE
)

if (pargs $ method == "cpm") {
  norm_factors <- 1000000 / colSums(expr_mat)
  expr_mat <- t(t(expr_mat) / norm_factors)

} else if (pargs $ method == "fpkm") {
  effleng <- pull(genes_meta, "merged")
  expr_mat <- apply(expr_mat, 2, function(counts) {
    exp(log(counts) + log(1e9) - log(effleng) - log(sum(counts)))
  })

} else if (pargs $ method == "tpm") {
  effleng <- pull(genes_meta, "merged")
  expr_mat <- apply(expr_mat, 2, function(counts) {
    rate <- log(counts) - log(effleng)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  })

} else if (pargs $ method == "uq") {
  norm_factors <- edgeR::calcNormFactors(
    object = expr_mat, lib.size = pull(sample_meta, "reads"),
    method = "upperquartile", p = pargs $ uq
  )
  expr_mat <- t(t(expr_mat) / norm_factors)

} else if (pargs $ method == "tmm") {
  norm_factors <- edgeR::calcNormFactors(
    object = expr_mat, lib.size = pull(sample_meta, "reads"),
    method = "tmm",
    logratioTrim = pargs $ m, sumTrim = pargs $ a,
    Acutoff = pargs $ ca, doWeighting = pargs $ weighting
  )
  expr_mat <- t(t(expr_mat) / norm_factors)

} else if (pargs $ method == "rle") {
  norm_factors <- edgeR::calcNormFactors(
    object = expr_mat, lib.size = pull(sample_meta, "reads"),
    method = "rle"
  )
  expr_mat <- t(t(expr_mat) / norm_factors)

} else if (pargs $ method == "seurat") {
  # do nothing.

} else if (pargs $ method == "sct") {
  # do nothing.

} else {
  cat("invalid method.", crlf)
  stop()
}

if (!dir.exists("norm"))
  dir.create("norm")

if (file.exists("norm/linear.rds"))
  file.remove("norm/linear.rds")
if (file.exists("norm/log.rds"))
  file.remove("norm/log.rds")

if (pargs $ method != "sct") {
  saveRDS(expr_mat, paste("norm/linear", pargs $ method, "rds", sep = "."))
  system(paste(
    "ln", "-s",
    paste(getwd(), paste("norm/linear", pargs $ method, "rds", sep = "."), sep = "/"), # nolint
    "norm/linear.rds"
  ))
}

# entering seurat.

dedup <- duplicated(pull(genes_meta, "seurat_names")) |
  pull(genes_meta, "seurat_names") == ""
expr_mat <- expr_mat[!dedup, ]
genes_meta <- genes_meta[!dedup, ]

rownames(expr_mat) <- pull(genes_meta, "seurat_names")
srat <- Seurat::CreateSeuratObject(
  counts = as(expr_mat, "sparseMatrix")
)

if (pargs $ method == "seurat" || pargs $ method == "sct") {
  srat <- Seurat::NormalizeData(
    srat,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
} else {
  logdata <- log1p(expr_mat)
  srat <- Seurat::SetAssayData(
    object = srat, layer = "data", new.data = logdata
  )
}

# cell cycle scoring:
#
# if we run a pca on our object, using the variable genes we found in
# findvariablefeatures() above, we see that while most of the variance can be
# explained by lineage, pc8 and pc10 are split on cell-cycle genes including
# top2a and mki67. we will attempt to regress this signal from the data,
# so that cell-cycle heterogeneity does not contribute to pca or
# downstream analysis.

# for methods other than sct, here, the srat object will contain two layers:
# counts and data already. we can now run the cell cycle scoring.

if (pargs $ cc) {
  s_genes <- taxo $ s
  g2m_genes <- taxo $ g2m
}

# first, we assign each cell a score, based on its expression of g2/m and s
# phase markers. these marker sets should be anticorrelated in their
# expression levels, and cells expressing neither are likely not cycling and
# in g1 phase.

# we assign scores in the cellcyclescoring() function, which stores s and
# g2/m scores in object meta data, along with the predicted classification
# of each cell in either g2m, s or g1 phase. cellcyclescoring() can also set
# the identity of the seurat object to the cell-cycle phase by passing
# set.ident = true (the original identities are stored as old.ident). please
# note that seurat does not use the discrete classifications (g2m/g1/s) in
# downstream cell cycle regression. instead, it uses the quantitative scores
# for g2m and s phase. however, we provide our predicted classifications in
# case they are of interest.

if (pargs $ cc) {
  srat <- Seurat::CellCycleScoring(
    srat, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE
  )

  sample_meta $ s_score <- srat @ meta.data $ S.Score
  sample_meta $ g2m_score <- srat @ meta.data $ G2M.Score
  sample_meta $ ccycle <- as.factor(srat @ meta.data $ Phase)
}


# set these named metadata to the seurat object

#. srat[["reads"]] <- pull(sample_meta, "reads")
#. srat[["detection"]] <- pull(sample_meta, "detection")
#. srat[["pct_ribo"]] <- pull(sample_meta, "pct_ribo")
#. srat[["pct_mito"]] <- pull(sample_meta, "pct_mito")
#. srat[["doublet_scores"]] <- pull(sample_meta, "doublet_scores")
#.
#. if (pargs $ cc) {
#.   srat[["s_score"]] <- pull(sample_meta, "s_score")
#.   srat[["g2m_score"]] <- pull(sample_meta, "g2m_score")
#. }

# we should set all the metadata columns, for later use in integration.

for (x in colnames(sample_meta)) {
  srat[[x]] <- pull(sample_meta, x)
}

if (pargs $ method == "sct") {
  srat <- Seurat::SCTransform(srat, vars.to.regress = pargs $ regress,
                              do.scale = TRUE, do.center = TRUE)
  sample_meta $ sct_reads <- srat @ meta.data $ nCount_SCT
  sample_meta $ sct_features <- srat @ meta.data $ nFeature_SCT
} else {
  srat <- Seurat::ScaleData(srat, vars.to.regress = pargs $ regress)
}

# file operations

if (pargs $ method == "seurat" || pargs $ method == "sct") {

  DefaultAssay(srat) <- "RNA"
  lognorm <- SeuratObject::LayerData(srat, layer = "data") |> as("sparseMatrix")
  saveRDS(lognorm, "norm/log.seurat.rds")

  if (pargs $ method == "seurat") {
    system(paste(
      "ln", "-s",
      paste(getwd(), paste("norm/log", pargs $ method, "rds", sep = "."), sep = "/"), # nolint
      "norm/log.rds"
    ))

  } else {

    DefaultAssay(srat) <- "SCT"
    lognorm <- SeuratObject::LayerData(srat, layer = "data") |>
      as("sparseMatrix")
    saveRDS(lognorm, "norm/log.sct.rds")
    corr <- SeuratObject::LayerData(srat, layer = "counts") |>
      as("sparseMatrix")
    saveRDS(corr, "norm/linear.sct.rds")

    #. ord <- c()
    #. for (cx in rownames(srat)) {
    #.   id <- which(genes_meta $ name == cx)
    #.   ord <- c(ord, id[1])
    #. }
    #.
    #. genes_meta <- genes_meta[ord, ]

    # seurat sc-transform will often (nearly always) wipe out some gene from
    # the sct dataset. however, we decide not to save the scaled.data object
    # the original rna dataset will not change.

    system(paste("ln", "-s", paste(getwd(), "norm/log.sct.rds", sep = "/"), "norm/log.rds")) # nolint
    system(paste("ln", "-s", paste(getwd(), "norm/linear.sct.rds", sep = "/"), "norm/linear.rds")) # nolint
  }

} else {
  lognorm <- SeuratObject::LayerData(srat, layer = "data") |> as("sparseMatrix")
  saveRDS(lognorm, paste("norm/log", pargs $ method, "rds", sep = "."))
  system(paste(
    "ln", "-s",
    paste(getwd(), paste("norm/log", pargs $ method, "rds", sep = "."), sep = "/"), # nolint
    "norm/log.rds"
  ))
}

genes_meta <- genes_meta |> tibble()

saveRDS(sample_meta, "norm/samples-meta.rds")
saveRDS(genes_meta, "norm/genes-meta.rds")
saveRDS(srat, paste("norm/seurat", pargs $ method, "rds", sep = "."))
system(paste(
  "ln", "-s",
  paste(getwd(), paste("norm/seurat", pargs $ method, "rds", sep = "."), sep = "/"), # nolint
  "norm/seurat.rds"
))

shared[["meta_sample"]] <- sample_meta
shared[["meta_gene"]] <- genes_meta
shared[["seurat"]] <- srat
