
require(reticulate)

scrublet <- function(
  seurat_obj,
  python_home               = Sys.which("python"),
  return_results_only       = FALSE,
  min_counts                = 2,
  min_cells                 = 3,
  expected_doublet_rate     = 0.06,
  min_gene_variability_pctl = 85,
  n_prin_comps              = 50,
  sim_doublet_ratio         = 2,
  n_neighbors               = NULL
) {

  reticulate::use_python(python_home)

  # test whether have the modules installed
  if (!reticulate::py_module_available("scrublet")) {
    stop("python module scrublet does not seem to be installed")
  }

  # source .py file
  reticulate::source_python(paste(gp_base, "scrublet.py", sep = "/"))

  # prepare the data
  x <- as(Matrix::t(seurat_obj @ assays $ RNA $ counts), "TsparseMatrix")
  i <- as.integer(x @ i)
  j <- as.integer(x @ j)
  val <- x @ x
  dim <- as.integer(x @ Dim)

  if (is.null(n_neighbors)) {
    n_neighbors <- round(0.5 * sqrt(nrow(x)))
  }

  # do the scrublet analysis
  scrublet_py_args <- c(list(
    i = i,
    j = j,
    val = val,
    dim = dim,
    expected_doublet_rate = expected_doublet_rate,
    min_counts = min_counts,
    min_cells = min_cells,
    min_gene_variability_pctl = min_gene_variability_pctl,
    n_prin_comps = n_prin_comps,
    sim_doublet_ratio = sim_doublet_ratio,
    n_neighbors = n_neighbors
  ))

  scrublet_res <- do.call(scrublet_py, scrublet_py_args) # nolint
  names(scrublet_res) <- c("doublet_scores", "predicted_doublets")

  if (return_results_only) {
    return(scrublet_res)

  } else {
    seurat_obj[["doublet_scores"]] <- scrublet_res $ doublet_scores
    seurat_obj[["predicted_doublets"]] <- scrublet_res $ predicted_doublets
    return(seurat_obj)
  }
}

scrublet_matrix <- function(
  seurat_obj, matrix,
  python_home               = Sys.which("python"),
  return_results_only       = FALSE,
  min_counts                = 2,
  min_cells                 = 3,
  expected_doublet_rate     = 0.06,
  min_gene_variability_pctl = 85,
  n_prin_comps              = 50,
  sim_doublet_ratio         = 2,
  n_neighbors               = NULL
) {

  reticulate::use_python(python_home)

  # test whether have the modules installed
  if (!reticulate::py_module_available("scrublet")) {
    stop("python module scrublet does not seem to be installed")
  }

  # source .py file
  reticulate::source_python(paste(gp_base, "scrublet.py", sep = "/"))

  # prepare the data
  x <- as(Matrix::t(matrix), "TsparseMatrix")
  i <- as.integer(x @ i)
  j <- as.integer(x @ j)
  val <- x @ x
  dim <- as.integer(x @ Dim)

  if (is.null(n_neighbors)) {
    n_neighbors <- round(0.5 * sqrt(nrow(x)))
  }

  # do the scrublet analysis
  scrublet_py_args <- c(list(
    i = i,
    j = j,
    val = val,
    dim = dim,
    expected_doublet_rate = expected_doublet_rate,
    min_counts = min_counts,
    min_cells = min_cells,
    min_gene_variability_pctl = min_gene_variability_pctl,
    n_prin_comps = n_prin_comps,
    sim_doublet_ratio = sim_doublet_ratio,
    n_neighbors = n_neighbors
  ))

  scrublet_res <- do.call(scrublet_py, scrublet_py_args) # nolint
  names(scrublet_res) <- c("doublet_scores", "predicted_doublets")

  if (return_results_only) {
    return(scrublet_res)

  } else {
    seurat_obj[["doublet_scores"]] <- scrublet_res $ doublet_scores
    seurat_obj[["predicted_doublets"]] <- scrublet_res $ predicted_doublets
    return(seurat_obj)
  }
}
