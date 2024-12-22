
require(ggplot2)
require(tidyr)

# the entropy of the genes across cell j is defined to be:
#
#     s1 = ln(mean(expr_j))
#
# this assumes that all cells in the cluster is homogenous. while the alternative
# assumption is that cells are totally heterogenous, and each cell represent its
# distinct state. the entropy for a totally heterogenous dataset is:
# 
#     s2 = mean(log(expr_ij))
#
# for calculation, this is simply a conversion of mean and log transform.
# r is a small positive value to avoid nan's been introduced due to invalid log.
# we require that the expression matrix have gene names as row names.

rogue.entropy <- function(expr, r = 1) {
  
  entropy <- Matrix::rowMeans(log(expr + 1))
  mean.expr <- log(Matrix::rowMeans(expr) + r)
  
  result <- tibble(
    gene = rownames(expr),
    mean.expr = mean.expr,
    entropy = entropy
  )
  
  return(result)
}


# fit the relationship between expression entropy and mean gene expression:
#
# in real-world dataset, there must be to some extent, heterogenous components
# in the dataset clusters. s1 > s2. and for different techniques for rna sequencing
# there is a systematic bias on the data distribution. we need a fitting step
# to make sure that the highly-variable genes in the clusters are picked.
#
#' @param span is a parameter for loess for controls of fit smoothness.
#' @param mt.method multiple testing method used in p.adjust.
#
# the returning dataframe consists of [1] gene names, [2] log mean expression (s1),
# [3] fit, the fitted curve of real world log mean expression trajectory (s1'),
# [4] ds, delta of entropy, the difference between s2 and s1',
# [5] p.value, assume that ds is normal, and the p value of normal distribution.
# [6] p.adj, the multi-test adjusted p value for normal approximation. and
# [7] entropy, the expression entropy per se.
#
# the returning values are sorted according to decreasing delta of entropy, so
# that the most highly-variable genes ranked on the top.

rogue.entropy.fit <- function(.x, span = 0.5, mt.method = "fdr"){
  
  .x <- .x %>%
    dplyr::filter(is.finite(mean.expr)) %>% 
    dplyr::filter(entropy > 0)
  
  fit <- loess(entropy ~ mean.expr, data = .x, span = span)
  prd <- predict(fit, .x $ mean.expr)
  .x %>%
    dplyr::mutate(fit = prd) %>%
    dplyr::mutate(ds = fit - entropy) %>%
    dplyr::mutate(pv = 1 - pnorm(. $ ds, mean = mean(. $ ds), sd = sd(. $ ds))) %>%
    dplyr::filter(pv > 0.1) -> tmp
  
  fit <- loess(entropy ~ mean.expr, data = tmp, span = span)
  prd <- predict(fit, .x $ mean.expr)
  .x %>%
    dplyr::mutate(fit = prd) %>%
    dplyr::mutate(ds = fit - entropy) %>%
    dplyr::filter(is.finite(ds)) %>%
    dplyr::mutate(pv = 1-pnorm(. $ ds, mean = mean(. $ ds), sd = sd(. $ ds))) %>%
    dplyr::filter(pv > 0.1) -> tmp
  
  fit <- loess(entropy ~ mean.expr, data = tmp, span = span)
  prd <- predict(fit, .x $ mean.expr)
  
  .x %>%
    dplyr::mutate(fit = prd) %>%
    dplyr::mutate(ds = fit - entropy) %>%
    dplyr::filter(is.finite(ds)) -> .x
  
  .x <- .x %>% dplyr::mutate(
    p.value = 1 - pnorm(.x $ ds, mean = mean(.x $ ds), sd = sd(.x $ ds)))
  
  p.adj <- p.adjust(.x $ p.value, method = mt.method)
  .x <- .x %>% dplyr::mutate(p.adj = p.adj) %>% dplyr::arrange(desc(ds))
}


# identify highly informative genes using rogue's s-e model
# 
# returns a tibble object with seven columns:
# [1] gene, the gene name.
# [2] mean.expr, the mean expression levels of genes.
# [3] entropy, the expected expression entropy from a given mean gene expression.
# [4] fit, the mean expression levels of genes.
# [5] ds, the entropy reduction against the null expectation.
# [6] p.value, the significance of ds against normal distribution.
# [7] p.adj, adjusted p value, or a copy of p.value if specifying not to adjust it.

rogue.hvg <- function(
    expr, span = 0.5, r = 1, 
    mt.method = "fdr", adjust.p.value = TRUE
  ) {
  
  entr <- rogue.entropy(expr, r = r)
  entr <- rogue.entropy.fit(entr, span = span, mt.method = mt.method)
  
  if(! adjust.p.value) { entr <- entr %>% dplyr::mutate(p.adj = p.value) }
  return(entr)
}


# filtering out low-abundance genes and low-quality cells

gene.cell.filter <- function(expr, min.cells = 10, min.genes = 10) {
  gene_count <- colSums(expr > 0, na.rm = T)
  cell_count <- rowSums(expr > 0, na.rm = T)
  lq1 <- cell_count < min.cells
  lq2 <- gene_count < min.genes
  return(expr[!lq1, !lq2])
}


# assess the purity of single cell population by rogue index.
#
# we may either specify the platform to access the default k value, or specify
# a k value manually. the rogue test is meaned to be run with highly variable
# genes. you can either specify them manually in `features`, or it is recommended
# to use the s-e model to pick the hvgs for you.

rogue.calc <- function(.x, platform = "umi", cutoff = 0.05, k = NULL, features = NULL) {
  
  if (is.null(k)) {
    if (is.null(platform)) {
      warning(
        "you should supply a constant k yourself, or specify the platform for ",
        "us to automatically pick the k constant for you. possible platforms ",
        "are `umi` or `full-length`, the default k value for each case is 45 ",
        "and 500 respectively."
      )
    } else if (platform == "umi"){
      k = 45
    } else if (platform == "full-length"){
      k = 500
    } else if (!is.null(platform) & !(platform %in% c("umi", "full-length"))) {
      warning("possible platforms: `umi` or `full-length`")
    }
    
  } else if (! is.null(k)) {
    k <- k
  }
  
  if (! is.null(features)) {
    .x <- .x %>% dplyr::filter(gene %in% features)
    sig_value <- sum(abs(.x $ ds))
    rogue <- 1 - sig_value / (sig_value + k)
    return(rogue)
    
  } else {
    sig_value <- abs(.x $ ds[.x $ p.adj < cutoff & .x $ p.value < cutoff])
    sig_value <- sum(sig_value)
    rogue <- 1 - sig_value / (sig_value + k)
    return(rogue)
  }
}


# remove outlier cells when calculating ROGUE
#' @param n remove this many outlier cells.

rogue.remove.outlier <- function(
    ent, expr, n = 2, span = 0.5, r = 1, mt.method = "fdr"
  ) {
  
  sig.gene <- ent %>% dplyr::filter(p.adj < 0.05) %>% dplyr::pull(gene)
  ng <- length(sig.gene)
  expr <- expr[sig.gene, ]
  
  mean.v <- c()
  entr.v <- c()
  for (i in 1:ng) {
    .x <- as.numeric(expr[i,])
    .x <- base::sort(.x, decreasing = T)
    .x <- .x[-c(1:n)]
    mean.v[i] <- log(mean(.x) + r)
    entr.v[i] <- mean(log(.x + 1))
  }
  
  mean.cut <- min(ent $ mean.expr)
  
  ent $ mean.expr[1:ng] <- mean.v
  ent $ entropy[1:ng] <- entr.v
  
  ent <- ent %>% dplyr::select(-p.adj) %>% dplyr::filter(mean.expr > mean.cut)
  ent <- rogue.entropy.fit(ent, span = span, mt.method = "fdr")
  return(ent)
}


# calculate rogue for each labels and samples in a dataset.
#
#' @param expr expression matrix. Rows should be genes and columns should be cells.
#' @param labels a vector of cell cluster labels
#' @param samples a vector of samples (e.g. patients) to which each cell belongs
#' @param min.cell.per.cluster only clusters containing at least this many cells will receive rogue values
#' @param remove.outlier.n remove this many out-lier cells when calculating rogue index
#' @param span parameter controls the degree of smoothing when fitting loess regression
#' @param r fixed value to avoid log(0) of mean gene expression levels.
#' @param filter whether to filter out low-abundance genes and low-quality cells.
#' @param mt.method multiple testing method used in adjusted p value.
#
# returns a matrix of rogue index with samples as rows, clusters as columns.

rogue <- function(
    expr, labels, samples, 
    platform = "umi", k = NULL,
    min.cell.per.cluster = 10, remove.outlier.n = 2, span = 0.5, r = 1,
    filter = FALSE, min.cells = 10, min.genes = 10,
    mt.method = "fdr"
  ) {
  
  clusters <- unique(labels)
  meta <- tibble(cell.id = 1 : ncol(expr), ct = labels, sample = samples)
  
  sample.rogue <- function(meta, cluster){
    tmp <- meta %>% dplyr::filter(ct == cluster)
    s <- unique(samples)
    rogue <- c()
    
    for (i in 1:length(s)) {
      index1 <- tmp %>% dplyr::filter(sample == s[i]) %>% dplyr::pull(cell.id)
      if (length(index1) >= min.cell.per.cluster) {
        tmp.matr <- expr[,index1]
        if (filter){
          print("filtering out low-abundance genes and low-quality cells ...")
          tmp.matr <- matr.filter(tmp.matr, min.cells = min.cells, min.genes = min.genes)
        } else { tmp.matr <- tmp.matr }
        
        tmp.res <- rogue.hvg(tmp.matr, span = span, r = r)
        tmp.res <- rogue.remove.outlier(tmp.res, tmp.matr, span = span, r = r, n = remove.outlier.n)
        rogue[i] <- rogue.calc(tmp.res, platform = platform, k = k)
      } else { rogue[i] <- NA }
    }
    
    return(rogue)
  }
  
  res <- list()
  
  for (i in 1:length(clusters)) {
    res[[i]] <- sample.rogue(meta, clusters[i])
  }
  
  res.tibble <- Reduce(rbind, res) %>% as.matrix() %>% t() %>% as.data.frame()
  colnames(res.tibble) <- clusters
  rownames(res.tibble) <- unique(samples)
  return(res.tibble)
}


# visualize the sample * cluster rogue matrix.

rogue.boxplot <- function(res.rogue){
  
  res.rogue %>%
    tidyr::gather(key = clusters, value = rogue.values) %>%
    ggplot(aes(clusters, rogue.values)) +
    geom_boxplot(color = "#FF3E96", outlier.shape = NA) +
    geom_point(color = "#FF3E96", size = 1.5) +
    theme_bw() +
    theme(axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 13, colour = "black")) +
    labs(
      x = "Clusters",
      y = "Rogue index"
    ) -> p
  
  return(p)
}


# calculate the value of the reference factor k

rogue.determine.k <- function(expr, span = 0.5, r = 1, mt.method = "fdr", adjust.p.value = TRUE) {
  
  res <- rogue.entropy(expr, r = r)
  res <- rogue.entropy.fit(res, span = span, mt.method = mt.method)
  
  if (! adjust.p.value) {
    res <- res %>% dplyr::mutate(p.adj = p.value)
  }
  
  k <- res %>% dplyr::filter(p.adj < 0.05) %>% dplyr::pull(ds) %>% sum()
  k <- k / 2
  return(k)
}
