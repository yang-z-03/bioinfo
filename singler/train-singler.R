
train_singler <- function(
    ref, labels,
    genes = "de",
    sd.thresh = NULL,
    de.method = c("classic", "wilcox", "t"),
    de.n = NULL,
    de.args = list(),
    aggr.ref = FALSE,
    aggr.args = list(),
    recompute = TRUE,
    restrict = NULL,
    assay.type = "logcounts",
    check.missing = TRUE,
    approximate = FALSE,
    num.threads = BiocParallel::bpnworkers(BPPARAM),
    BNPARAM = NULL,
    BPPARAM = BiocParallel::SerialParam()) {

  de.method <- match.arg(de.method)

  if (solo <- !.is_list(ref)) {
    ref <- list(ref)
    labels <- list(labels)
    if (!is.character(genes)) {
      genes <- list(genes)
    }
  }

  if (!BiocParallel::bpisup(BPPARAM) && !is(BPPARAM, "MulticoreParam")) {
    BiocParallel::bpstart(BPPARAM)
    on.exit(BiocParallel::bpstop(BPPARAM))
  }

  ref <- lapply(ref,
    FUN = .to_clean_matrix, assay.type = assay.type,
    check.missing = check.missing, msg = "ref", BPPARAM = BPPARAM
  )

  # cleaning the genes.
  gns <- lapply(ref, rownames)
  if (length(unique(gns)) != 1L) {
    stop("row names are not identical across references")
  }

  if (!is.null(restrict)) {
    keep <- gns[[1]] %in% restrict
    ref <- lapply(ref, FUN = "[", i = keep, , drop = FALSE)
  }

  if (S4Vectors::isSingleString(genes)) {
    genes <- rep(genes, length(ref))
  } else if (length(genes) != length(ref)) {
    stop("list-like 'genes' should be the same length as 'ref'")
  }

  # cleaning the labels.
  labels <- lapply(labels, as.character)
  if (length(labels) != length(ref)) {
    stop("lists in 'labels' and 'ref' should be of the same length")
  }

  for (l in seq_along(labels)) {
    keep <- !is.na(labels[[l]])
    if (!all(keep)) {
      labels[[l]] <- labels[[l]][keep]
      ref[[l]] <- ref[[l]][, keep, drop = FALSE]
    }
  }

  gene.info <- mapply(
    FUN = .identify_genes, ref = ref, labels = labels, genes = genes,
    MoreArgs = list(
      de.method = de.method, de.n = de.n,
      de.args = de.args, BPPARAM = BPPARAM),
    SIMPLIFY = FALSE
  )

  if (!solo && !recompute) {
    .Deprecated(old = "recompute = FALSE")
  }

  output <- mapply(
    FUN = .build_trained_index, ref = ref, labels = labels, markers = gene.info,
    MoreArgs = list(
      aggr.ref = aggr.ref, aggr.args = aggr.args, BPPARAM = BPPARAM,
      approximate = approximate, num.threads = num.threads),
    SIMPLIFY = FALSE
  )

  if (solo) {
    output[[1]]
  } else {
    final <- S4Vectors::List(output)
    S4Vectors::metadata(final) $ recompute <- recompute
    final
  }
}

.identify_genes <- function(
    ref, labels, genes = "de", de.method = "classic", de.n = NULL,
    de.args = list(), BPPARAM = BPPARAM) {

  if (length(labels) != ncol(ref)) {
    stop("number of labels must be equal to number of cells")
  }

  if (.is_list(genes)) {
    is.char <- vapply(genes, is.character, TRUE)
    if (all(is.char)) {
      genes <- .convert_per_label_set(genes)
    } else if (any(is.char)) {
      stop("'genes' must be a list of character vectors or a list of list of vectors")
    }

    # to convert from S4Vectors::List of S4Vectors::Lists.
    genes <- lapply(genes, as.list) 
    genes <- .validate_de_gene_set(genes, labels)

    # ensure that the user hasn't supplied genes that aren't available.
    rn <- rownames(ref)
    genes <- lapply(genes, function(l) lapply(l, intersect, rn))

  } else {
    genes <- match.arg(genes, c("de", "sd", "all"))
    if (genes != "de") {
      .Deprecated(old = "genes = \"", genes, "\"")
    }
    genes <- .get_genes_by_de(
      ref, labels, de.n = de.n, de.method = de.method,
      de.args = de.args, BPPARAM = BPPARAM
    )
  }

  genes
}

.build_trained_index <- function(
    ref, labels, markers, aggr.ref, aggr.args, search.info, 
    approximate = FALSE, BPPARAM = BiocParallel::SerialParam(), num.threads = 1) {

  if (aggr.ref) {
    aggr <- do.call(aggregate_reference, c(list(
      ref = quote(ref), label = labels,
      check.missing = FALSE, BPPARAM = BPPARAM
    ), aggr.args))
    ref <- assay(aggr)
    labels <- aggr $ label
  }

  if (anyNA(labels)) {
    keep <- !is.na(labels)
    ref <- ref[, keep, drop = FALSE]
    labels <- labels[keep]
  }

  ulabels <- .get_levels(labels)

  built <- .build_index(
    ref, markers = markers, labels = labels, ulabels = ulabels,
    approximate = approximate, num.threads = num.threads
  )

  S4Vectors::List(
    built = built, ref = ref,
    labels = list(full = labels, unique = ulabels),
    markers = list(full = markers, unique = rownames(ref)[get_subset(built) + 1]),
    options = list(approximate = approximate)
  )
}

.build_index <- function(
    ref, markers, labels, ulabels, approximate, num.threads ) {

  for (m in seq_along(markers)) {
    current <- markers[[m]]
    for (n in seq_along(current)) {
      idx <- match(current[[n]], rownames(ref))
      if (anyNA(idx)) {
        stop("could not find '", 
             current[[n]][which(is.na(idx))[1]], 
             "' in 'rownames(ref)'")
      }
      current[[n]] <- idx - 1L
    }
    markers[[m]] <- current
  }

  parsed <- beachmat::initializeCpp(ref)
  prebuild(parsed, match(labels, ulabels) - 1L, markers,
           approximate = approximate, nthreads = num.threads)
}

.get_levels <- function(labels) sort(unique(labels))

# Unfortunately, we can't test for List, because each trained structure is
# also a list; so we just check whether the 'ref' field exists.
.is_solo <- function(trained) !is.null(trained $ ref)

.get_genes_by_de <- function(
    ref, labels,
    de.method = "classic", de.n = NULL, de.args = list(),
    BPPARAM = BiocParallel::SerialParam()
  ) {

  if (de.method == "classic") {
    get_classic_markers(ref = ref, labels = labels, 
                      de.n = de.n, check.missing = FALSE,
                      BPPARAM = BPPARAM)
  } else {
    if (de.method == "t") {
      FUN <- scran::pairwiseTTests
    } else {
      FUN <- scran::pairwiseWilcox
    }

    pairwise <- do.call(FUN, c(list(
      x = ref, groups = labels, direction = "up",
      log.p = TRUE, BPPARAM = BPPARAM), de.args))
    
    if (is.null(de.n)) {
      de.n <- 10
    }

    collected <- scran::getTopMarkers(
      pairwise $ statistics, pairwise $ pairs, n = de.n,
      pval.field = "log.p.value", fdr.field = "log.FDR",
      fdr.threshold = log(0.05)
    )

    lapply(collected, as.list)
  }
}

.convert_per_label_set <- function(genes) {
  
  # Converting into a list-of-lists format so that it plays nice with
  # downstream methods. This is done by saving the markers on the diagonal
  # so that each label's markers are included in the set of genes to use 
  # during fine-tuning. Don't exclude the diagonal!

  all.labs <- names(genes)
  for (i in all.labs) {
    empty <- rep(list(character(0)), length(all.labs))
    names(empty) <- all.labs
    empty[[i]] <- genes[[i]]
    genes[[i]] <- empty
  }
  genes
}

.validate_de_gene_set <- function(genes, labels) {
  ulabels <- .get_levels(labels)
  if (!all(ulabels %in% names(genes))) {
    stop("need marker gene information for each label")
  }

  genes <- genes[ulabels]
  for (u in ulabels) {
    if (!all(ulabels %in% names(genes[[u]]))) {
      stop("need marker genes between each pair of labels")
    }
    genes[[u]] <- genes[[u]][ulabels]
  }

  genes
}
