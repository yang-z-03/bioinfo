
require(monocle)

# 1.  LOADING THE DATASET ======================================================

## 1.1  the universal solution -------------------------------------------------

# Monocle holds single cell expression data in objects of the CellDataSet class.
# The class is derived from the Bioconductor ExpressionSet class, which provides
# a common interface familiar to those who have analyzed microarray experiments
# with Bioconductor. The class requires three input files:

# - exprs, a numeric matrix of expression values, where rows are genes,
#   and columns are cells
# - phenoData, an AnnotatedDataFrame object, where rows are cells, and columns
#   are cell attributes (such as cell type, culture condition, day captured)
# - featureData, an AnnotatedDataFrame object, where rows are features
#   (e.g. genes), and columns are gene attributes, such as biotype, gc content

# Required dimensions for input files
# The expression value matrix must:
# - have the same number of columns as the phenoData has rows.
# - have the same number of rows as the featureData data frame has rows.
#
# Additionally:
# - row names of the phenoData object should match the column names of
#   the expression matrix.
# - row names of the featureData object should match row names of the
#   expression matrix.
# - one of the columns of the featureData should be named "gene_short_name".

# You should import only the RAW counts (UMI counts if your experiment is UMI-
# based, or FKPM or TPM if not). Do not normalize the data by yourself, since
# monocle will do this for you.

pd <- new("AnnotatedDataFrame", data = shared[["meta_sample"]])
fd <- new("AnnotatedDataFrame", data = shared[["meta_gene"]])

# Monocle works well with both relative expression data and count-based measures
# (e.g. UMIs). In general, it works best with transcript count data, especially
# UMI data. Whatever your data type, it is critical that specify the appropriate
# distribution for it. FPKM/TPM values are generally log-normally distributed,
# while UMIs or read counts are better modeled with the negative binomial. To
# work with count data, specify the negative binomial distribution as the
# `expressionFamily` argument.

# family function       note
#
# negbinomial.size()    UMIs, Transcript counts from experiments with spike-ins
#                       or relative2abs(), raw read counts
#                       Negative binomial distribution with fixed variance
#                       (which is automatically calculated by Monocle).
#                       Recommended for most users.
#
# negbinomial()         UMIs, Transcript counts from experiments with spike-ins
#                       or relative2abs, raw read counts
#                       Slightly more accurate than negbinomial.size(), but much
#                       much slower. Not recommended except for very small
#                       datasets.
#
# tobit()               FPKM, TPM
#                       Tobits are truncated normal distributions. Using tobit()
#                       will tell Monocle to log-transform your data where
#                       appropriate. Do not transform it yourself.
#
# gaussianff()          log-transformed FPKM/TPMs, Ct values from single-cell
#                       qPCR If you want to use Monocle on data you have already
#                       transformed to be normally distributed, you can use
#                       this function, though some Monocle features may
#                       not work well.

# avoid the use of as.matrix(), because this turns your matrix into a dense one.
# and will take up a great number of memory!

cds <- newCellDataSet(
  as(counts(shared[["seurat"]]), "sparseMatrix"),
  phenoData = pd, featureData = fd
)


## 1.2  read from 10x cell ranger ----------------------------------------------

# If you have 10X Genomics data and are using cellrangerRkit, you can use it
# to load your data and then pass that to Monocle as follows:

cellranger_pipestance_path <- getwd()
gbm <- load_cellranger_matrix(cellranger_pipestance_path)

fd <- fData(gbm)

# The number 2 is picked arbitrarily in the line below.
# Where "2" is placed you should place the column number that corresponds
# to your featureData's gene short names.

colnames(fd)[2] <- "gene_short_name"

gbm_cds <- newCellDataSet(
  exprs(gbm), # already a sparse matrix.
  phenoData = new("AnnotatedDataFrame", data = pData(gbm)),
  featureData = new("AnnotatedDataFrame", data = fd),
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size()
)


## 1.3  converting to absolute expressions
##      with or without a spike-in test    -------------------------------------

# If you performed your single-cell RNA-Seq experiment using spike-in standards,
# you can convert these measurements into mRNAs per cell (RPC). RPC values are
# often easier to analyze than FPKM or TPM values, because have better
# statistical tools to model them. In fact, it's possible to convert FPKM or TPM
# values to RPC values even if there were no spike-in standards included in the
# experiment. Monocle 2 includes an algorithm called Census which performs this
# conversion. You can convert to RPC values before creating your CellDataSet
# object using the relative2abs() function, as follows:

pd <- new("AnnotatedDataFrame", data = shared[["meta_sample"]])
fd <- new("AnnotatedDataFrame", data = shared[["meta_gene"]])

# first create a celldataset from the relative expression levels
expr <- newCellDataSet(
  as.matrix(counts(shared[["seurat"]])), # is this matrix dense?
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.1,
  expressionFamily = tobit(Lower = 0.1)
)

# next, use it to estimate rna counts
rpc_matrix <- relative2abs(expr, method = "num_genes")

# now, make a new celldataset using the rna counts
expr <- newCellDataSet(
  as(as.matrix(rpc_matrix), "sparseMatrix"),
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.5, # note this alters!
  expressionFamily = negbinomial.size()
)


## 1.4  normalization ----------------------------------------------------------

# Finally, we'll also call two functions that pre-calculate some information
# about the data. Size factors help us normalize for differences in mRNA
# recovered across cells, and "dispersion" values will help us perform
# differential expression analysis later.

# [!] NOTE:
#
#   estimateSizeFactors() and estimateDispersions() will only work, and are
#   only needed, if you are working with a CellDataSet with a negbinomial()
#   or negbinomial.size() expression family.

expr <- estimateSizeFactors(expr)
expr <- estimateDispersions(expr)
