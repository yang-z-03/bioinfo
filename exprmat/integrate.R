
parser <- argparse::ArgumentParser(
  prog = "integrate",
  description = "integrate qc-ed and normalized subset datasets into one"
)

parser $ add_argument(
  "-i", dest = "input", type = "character", nargs = "*", default = c(),
  help = paste("input sample directory name")
)

parser $ add_argument(
  "-m", dest = "method", type = "character",
  help = paste("method of normalization. supports 'cca', 'rpca', 'jpca',",
               "'harmony', 'scvi', and 'mnn'.")
)

parser $ add_argument(
  "-d", dest = "dataset", type = "character",
  help = paste("the column indicating dataset identity (the batches we want",
               "to integrate) within the sample metadata.")
)

parser $ add_argument(
  "--orig", dest = "origin", type = "character", default = "pca",
  help = paste("original reduction, by default 'pca'")
)

parser $ add_argument(
  "--dim", dest = "dim", type = "integer", default = 30,
  help = paste("pca dimension components taken into calculation")
)

parser $ add_argument(
  "--int", dest = "integration", type = "character", default = "",
  help = paste("integration reduction name. by default it will be set to the",
               "method name prefixed with 'integrated.'")
)

parser $ add_argument(
  "--n-feat", dest = "nfeature", type = "integer", default = 3000,
  help = paste("number of features to used for anchor construction. larger",
               "values may result in more memory consumption and slower speed")
)

parser $ add_argument(
  "--add-feat", dest = "addition", type = "character", nargs = "*", default = c(), # nolint
  help = paste("additional feature names (gene names) passed as features,",
               "specified manually, but required to exist in all the samples")
)

parser $ add_argument(
  "--k-anchor", dest = "k.anchor", type = "integer", default = 10,
  help = paste("number of neighbors used when picking anchors")
)

parser $ add_argument(
  "--k-filter", dest = "k.filter", type = "integer", default = 200,
  help = paste("number of neighbors used when filtering anchors")
)

parser $ add_argument(
  "--merge-only", dest = "merge", action = "store_true", default = FALSE,
  help = "only merge the dataset, without any integration mapping"
)

parser $ add_argument(
  "--sct", dest = "sct", action = "store_true", default = FALSE,
  help = "run the integration with SCT normalization"
)

parser $ add_argument(
  "--use-sct-corrected", dest = "sctc", action = "store_true", default = FALSE,
  help = "read the SCT-corrected data for each sample rather than raw counts"
)

parser $ add_argument(
  "--no-norm", dest = "nonorm", action = "store_true", default = FALSE,
  help = "do not normalize"
)

parser $ add_argument(
  "--mnn-k", dest = "mnn.k", type = "integer", default = 20,
  help = paste("an integer scalar specifying the number of nearest neighbors",
               "to consider when identifying mnns")
)

parser $ add_argument(
  "--mnn-dim", dest = "mnn.dim", type = "integer",
  default = 50,
  help = paste("numeric scalar specifying the number of dimensions to use for",
               "dimensionality reduction in multiBatchPCA. If NA, no",
               "dimensionality reduction is performed and any matrices",
               "in ... are used as-is")
)

parser $ add_argument(
  "--mnn-dist", dest = "mnn.dist", type = "integer", default = 3,
  help = paste("a numeric scalar specifying the threshold beyond which",
               "neighbours are to be ignored when computing correction",
               "vectors. each threshold is defined as a multiple of the",
               "number of median distances")
)

parser $ add_argument(
  "--mnn-pk", dest = "mnn.pk", type = "double", default = -1,
  help = paste("a numeric scalar in (0, 1) specifying the proportion of cells",
               "in each dataset to use for mutual nearest neighbor searching.",
               "if set, the number of nearest neighbors used for the MNN",
               "search in each batch is redefined as max(k, pk * N) where N is",
               "the number of cells in that batch")
)

parser $ add_argument(
  "--scvi-dim", dest = "scvi.dim", type = "integer",
  default = 30,
  help = paste("dimensionality of the latent space")
)

parser $ add_argument(
  "--scvi-layer", dest = "scvi.layer", type = "integer", default = 2,
  help = paste("number of hidden layers used for encoder and decoder NNs.")
)

parser $ add_argument(
  "--scvi-gdist", dest = "scvi.gdist", type = "character", default = "nb",
  help = paste("one of: 'nb' - negative binomial distribution, 'zinb' -",
               "zero-inflated negative binomial distribution, 'poisson' -",
               "poisson distribution")
)

parser $ add_argument(
  "--scvi-epoch", dest = "scvi.epoch", type = "integer", default = -1,
  help = paste("maximum training epoch")
)

parser $ add_argument(
  "--regress", dest = "regr", type = "character", nargs = "*", default = c(),
  help = paste("variables to regress when scaling")
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

if (pargs $ mnn.pk < 0) pargs $ mnn.pk <- NA
if (pargs $ scvi.epoch < 0) pargs $ scvi.epoch <- NULL
pargs $ scvi.dim <- min(pargs $ scvi.dim, pargs $ dim)
pargs $ mnn.dim <- min(pargs $ mnn.dim, pargs $ dim)

# be aware that so far there is no 'best' integration method for all scenarios.
# it is therefore important to try different methods and compare, to at the end
# choose the one that works the best for every specific case.

# in previous versions of seurat, we would require the data to be represented as
# two different seurat objects. in seurat v5, we keep all the data in one
# object, but simply split it into multiple ‘layers’.

# the method

norm <- "LogNormalize"
if (pargs $ sct) norm <- "SCT"

seurat_raw_data_only <- function(srat, pargs) {
  # the meta.data items i want to keep
  idents_df <- srat @ meta.data
  assayname <- "RNA"
  if (pargs $ sctc) assayname <- "SCT"

  export <- Seurat::CreateSeuratObject(
    counts = SeuratObject::LayerData(
      srat, assay = assayname,
      layer = "counts"
    ),
    data = SeuratObject::LayerData(
      srat, assay = assayname,
      layer = "data"
    ),
    project = "integrated",
    meta.data = idents_df
  )

  return(export)
}

# we should first merge our datasets. and this generates features already
# splitted by runs. however, considering the need to manually specify grouping
# information, we will merge the default splits, and re-split by ourself
# make sure you have a shared meta data column indicating the grouping.

if (file.exists("norm/seurat.rds")) {
  merged <- readRDS("norm/seurat.rds")

} else {

  if (length(pargs $ input) > 0) {
    run_names <- pargs $ input
  } else {
    run_names <- list.dirs(".", full.names = FALSE, recursive = FALSE)
  }

  seurat_list <- list()
  for (runname in run_names) {
    if (runname == "data") next
    if (runname == "norm") next

    if (file.exists(paste(".", runname, "norm", "seurat.rds", sep = "/"))) {
      original <- readRDS(
        paste(".", runname, "norm", "seurat.rds", sep = "/")
      )

      # read the metadata into the seurat object.
      #. assayname <- "RNA"
      #. if (pargs $ sctc) assayname <- "SCT"
      #. genemeta_name <- "genes-meta.rds"
      #. if (pargs $ sctc) genemeta_name <- "sct-genes-meta.rds"
      #. geneinfo <- readRDS(
      #.   paste(".", runname, "norm", genemeta_name, sep = "/")
      #. )
      #.
      #. print(geneinfo)

      seurat_list[[runname]] <- original |> seurat_raw_data_only(pargs)

      #. for (gcol in colnames(geneinfo)) {
      #.   vals <- pull(geneinfo, gcol)
      #.   names(vals) <- pull(geneinfo, "name")
      #.   seurat_list[[runname]][[assayname]][[gcol]] <- pull(geneinfo, gcol)
      #. }

      rm(original)

    } else {
      cat(red(runname), "does not contain a norm/seurat.rds", crlf)
      cat(yellow("make sure you run the norm procedure for each sample.\n"))
      stop()
    }
  }

  rmv <- c("data", "norm")
  run_names <- run_names[!(run_names %in% rmv)]

  # simply merge the list

  if (length(run_names) == 0) {
    cat(yellow("no sample directories."), crlf)
    stop()
  }

  if (length(run_names) == 1) {
    merged <- seurat_list[[run_names[1]]]
  } else {
    merged <- seurat_list[[run_names[1]]]
    cat(yellow("merging"), run_names[1], crlf)
    seurat_list[[run_names[1]]] <- NULL
    merged <- merge(merged, seurat_list, add.cell.ids = TRUE)
  }

  rm(seurat_list)
  suppressMessages(gc())

  # join the automatic splits

  merged[["RNA"]] <- JoinLayers(merged[["RNA"]])
}

if (pargs $ merge) {
  saveRDS(merged, "norm/seurat.rds")
  stop()
}

# re-split myself

# you have scaled your data prior to the integration, so we will run PCA
# directly, this require you pick the variable features and scale the data

options(future.globals.maxSize = 1000 * 1024^2)
cat(blue("running pca for unintegrated data ..."), crlf)

merged[["RNA"]] <- split(merged[["RNA"]],
                         f = merged @ meta.data[[pargs $ dataset]])

if (!pargs $ nonorm) {
  if (!pargs $ sct) {
    merged <- ScaleData(
      merged, vars.to.regress = c(pargs $ dataset, pargs $ regr)
    )
    merged <- FindVariableFeatures(merged)
  } else {
    # here we do not use the variables to regress ...
    merged <- Seurat::SCTransform(merged, do.scale = TRUE, do.center = TRUE)
  }
}

merged <- RunPCA(merged, npcs = pargs $ dim, verbose = FALSE)
merged <- RunUMAP(merged, dims = 1 : pargs $ dim, reduction = "pca",
                  reduction.name = "umap.unintegrated")

# here comes the parameter specification and integration process.

methodfunc <- NULL

cat(crlf, blue("selecting integration features ..."), crlf)
features <- Seurat::SelectIntegrationFeatures(
  list(merged), nfeatures = pargs $ nfeature,
  fvf.nfeatures = pargs $ nfeature
)

features <- c(features, pargs $ addition)
features <- features[!duplicated(features)]

redname <- pargs $ integration
if (redname == "") redname <- paste("integrated", pargs $ method, sep = ".")

cat(crlf, blue("integrating layers ..."), crlf)

if (pargs $ method == "cca") {
  methodfunc <- CCAIntegration
  merged <- IntegrateLayers(
    object = merged,
    method = methodfunc,
    orig.reduction = pargs $ origin,
    new.reduction = redname,
    normalization.method = norm,
    dims = 1 : pargs $ dim,
    features = features,
    k.anchor = pargs $ k.anchor,
    k.filter = pargs $ k.filter,
    verbose = TRUE
  )

} else if (pargs $ method == "rpca") {
  methodfunc <- RPCAIntegration
  merged <- IntegrateLayers(
    object = merged,
    method = methodfunc,
    orig.reduction = pargs $ origin,
    new.reduction = redname,
    normalization.method = norm,
    dims = 1 : pargs $ dim,
    features = features,
    k.anchor = pargs $ k.anchor,
    k.filter = pargs $ k.filter,
    verbose = TRUE
  )

} else if (pargs $ method == "mnn") {
  methodfunc <- FastMNNIntegration
  merged <- IntegrateLayers(
    object = merged,
    method = methodfunc,
    orig.reduction = pargs $ origin,
    new.reduction = redname,
    normalization.method = norm,
    dims = 1 : pargs $ dim,
    features = features,
    k.anchor = pargs $ k.anchor,
    k.filter = pargs $ k.filter,
    verbose = TRUE,

    # arguments passed to fastmnn
    k = pargs $ mnn.k,
    prop.k = pargs $ mnn.pk,
    ndist = pargs $ mnn.dist,
    d = pargs $ mnn.dim
  )

} else if (pargs $ method == "harmony") {
  methodfunc <- HarmonyIntegration
  merged <- IntegrateLayers(
    object = merged,
    method = methodfunc,
    orig.reduction = pargs $ origin,
    new.reduction = redname,
    normalization.method = norm,
    dims = 1 : pargs $ dim,
    features = features,
    k.anchor = pargs $ k.anchor,
    k.filter = pargs $ k.filter,
    verbose = TRUE,

    # arguments passed to harmony
    npcs = pargs $ dim
  )

} else if (pargs $ method == "scvi") {
  methodfunc <- scVIIntegration
  merged <- IntegrateLayers(
    object = merged,
    method = methodfunc,
    orig.reduction = pargs $ origin,
    new.reduction = redname,
    normalization.method = norm,
    dims = 1 : pargs $ dim,
    features = features,
    k.anchor = pargs $ k.anchor,
    k.filter = pargs $ k.filter,
    verbose = TRUE,

    # arguments passed to harmony
    ndims = pargs $ scvi.dim,
    nlayers = pargs $ scvi.layer,
    gene_likelihood = pargs $ gdist,
    max_epochs = pargs $ epoch
  )

} else if (pargs $ method == "jpca") {
  methodfunc <- JointPCAIntegration
  merged <- IntegrateLayers(
    object = merged,
    method = methodfunc,
    orig.reduction = pargs $ origin,
    new.reduction = redname,
    normalization.method = norm,
    dims = 1 : pargs $ dim,
    features = features,
    k.anchor = pargs $ k.anchor,
    k.filter = pargs $ k.filter,
    verbose = TRUE
  )
}

merged[["RNA"]] <- JoinLayers(merged[["RNA"]])

cat(crlf, blue("running umap dimension reduction ..."), crlf)
umapname <- paste("umap", pargs $ method, sep = ".")
merged <- RunUMAP(merged, reduction = redname, dims = 1 : pargs $ dim,
                  reduction.name = umapname)

if (!dir.exists("norm"))
  dir.create("norm")

saveRDS(merged, "norm/seurat.rds")

shared[["seurat"]] <- merged

# The approach of first normalizing each sample (matrix) is only advisable if
# your samples have roughly the same celltype compositions and you want to
# remove batch effects that are characterized by simple shifts in mean
# expression.
#
# -- https://github.com/satijalab/sctransform/issues/55#issuecomment-633843730

# by default we will merge the raw count as well as the log transformed data.
# and ignore the already-scaled layer in the assay. and run scale again here.
# for sct assays, we need to specify an extra tag --use-sct-corrected to include
# sct-corrected counts and log-space data in. by default, we will read the
# raw umi counts from the RNA assay and its log counts (run with log normalize)
# in norm.R. and run a SCT in the whole merge here.

# running SCT separately for samples is not recommended since it may lose
# information of different cell clusters. you need to make sure that your cells
# and experimental designs are completely homologous before setting the
# --use-sct-corrected here. (see above github discussion)

# run intgmeta to generate a combined gene metadata file

if (length(pargs $ input) > 0) {
  run_names <- pargs $ input
} else {
  run_names <- list.dirs(".", full.names = FALSE, recursive = FALSE)
}

merged_gene_info <- NULL
for (runname in run_names) {
  if (runname == "data") next
  if (runname == "norm") next

  geneinfo <- readRDS(
    paste(".", runname, "norm", "genes-meta.rds", sep = "/")
  )

  if (merged_gene_info |> is.null()) {
    merged_gene_info <- geneinfo
  } else {
    merged_gene_info <- dplyr::bind_rows(merged_gene_info, geneinfo)
  }
}

dup <- duplicated(pull(merged_gene_info, "seurat_names")) |
  pull(merged_gene_info, "seurat_names") == ""

merged_gene_info <- merged_gene_info[!dup, ]

# filtering and ordering

ord <- c()
for (cx in rownames(shared[["seurat"]])) {
  id <- which(merged_gene_info $ seurat_names == cx)
  ord <- c(ord, id[1])
}

meta <- merged_gene_info[ord, ]
saveRDS(meta, "norm/genes-meta.rds")
shared[["meta_gene"]] <- meta
