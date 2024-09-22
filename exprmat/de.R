
parser <- argparse::ArgumentParser(
  prog = "de",
  description = "get a list of differential expression genes amoung clusters"
)

parser $ add_argument(
  "-m", dest = "method", type = "character", default = "wilcox",
  help = paste("method of calculating differences. one of 'wilcox', 'limma',",
               "'bimod', 'roc', 't', 'negbinom', 'poisson', 'logistic',",
               "'mast', or 'deseq2'")
)

parser $ add_argument(
  "--cmp", dest = "cmp.cluster",
  type = "character", default = c(), nargs = "+",
  help = paste("cluster names (combined with OR) as the comparing set")
)

parser $ add_argument(
  "--with", dest = "with.cluster",
  type = "character", default = c(),
  help = paste("cluster names (combined with OR) as the reference set")
)

parser $ add_argument(
  "--cmpc", dest = "cmp.cond",
  type = "character", default = c(), nargs = "+",
  help = paste("conditions specifying the comparing set")
)

parser $ add_argument(
  "--withc", dest = "with.cond",
  type = "character", default = c(), nargs = "+",
  help = paste("conditions as the reference set. note that --cmpc and --withc",
               "must be used simutaneously to replace the use of --cmp and",
               "--with. this has higher precedency so that specification of",
               "this pair of options will overwrite those cluster names",
               "specified by --cmp and --with. the grammar of condition",
               "specification is the same as those required in dimreduc")
)

parser $ add_argument(
  "--top", dest = "top", type = "integer", default = 10,
  help = paste("only return the top N de-genes")
)

parser $ add_argument(
  "--pseudo-bulk", dest = "psbulk", default = c(), nargs = "*",
  help = paste("make pseudo-bulk references by classifying through these",
               "specified categories. note that the conditions to compare",
               "clusters must set within the pseudo-bulking variables. we",
               "will not check this, and the program will simply fail or",
               "behave wierd if you do not follow this. we recommend using",
               "--cmpc and --withc, but you can use --cmp and --with also if",
               "you put the identity column into the list.")
)

parser $ add_argument(
  "--min-pct", dest = "mpct", type = "double", default = 0.01,
  help = paste("threshold for minimum marker precentage")
)

parser $ add_argument(
  "--min-pct-d", dest = "mdpct", type = "double", default = -1,
  help = paste("threshold for minimum marker precentage difference between",
               "comparing and compared set")
)

parser $ add_argument(
  "--max-c", dest = "maxc", type = "integer", default = 1000000,
  help = paste("downsample the clusters with too much cells")
)

parser $ add_argument(
  "--min-logfc", dest = "mfc", type = "double", default = 0.5,
  help = paste("threshold for minimum log fold expression change")
)

parser $ add_argument(
  "--p", dest = "p", type = "double", default = 0.05,
  help = paste("threshold for p-value screen, set to 1 for no screen")
)

parser $ add_argument(
  "--store", dest = "store", type = "character", default = "",
  help = paste("store the de gene list to the shared objects")
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

auto_cast <- function(binary) {
  destclass <- shared[[".psbulk"]] @ meta.data |>
    dplyr::pull(binary[1]) |>
    class()

  if (destclass == "factor")
    destclass <- shared[[".psbulk"]] @ meta.data |>
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
    return(shared[[".psbulk"]][[binary[1]]] >= comp)

  } else if (stringr::str_detect(cond, "<=")) {
    binary <- stringr::str_split_1(cond, "<=")
    comp <- auto_cast(binary)
    return(shared[[".psbulk"]][[binary[1]]] <= comp)

  } else if (stringr::str_detect(cond, ">")) {
    binary <- stringr::str_split_1(cond, ">")
    comp <- auto_cast(binary)
    return(shared[[".psbulk"]][[binary[1]]] > comp)

  } else if (stringr::str_detect(cond, "<")) {
    binary <- stringr::str_split_1(cond, "<")
    comp <- auto_cast(binary)
    return(shared[[".psbulk"]][[binary[1]]] < comp)

  } else if (stringr::str_detect(cond, "=")) {
    binary <- stringr::str_split_1(cond, "=")
    comp <- auto_cast(binary)
    return(shared[[".psbulk"]][[binary[1]]] == comp)
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

find_marker <- function(pargs, obj, ident.1, ident.2, verbose) {
  FindMarkers(
    obj, ident.1 = ident.1, ident.2 = ident.2, verbose = verbose,
    min.pct = pargs $ mpct, min.diff.pct = pargs $ mdpct,
    max.cells.per.ident = pargs $ maxc,
    logfc.threshold = pargs $ mfc,
    min.cells.group = 1,
    test.use = switch(
      pargs $ method,
      wilcox = "wilcox",
      limma = "wilcox_limma",
      bimod = "bimod",
      roc = "roc",
      t = "t",
      negbinom = "negbinom",
      poisson = "poisson",
      logistic = "LR",
      mast = "MAST",
      deseq2 = "DESeq2"
    )
  )
}

# generate the pseudo-bulk equility

if (pargs $ psbulk |> length() > 0) {
  shared[[".psbulk"]] <- AggregateExpression(
    shared[["seurat"]], assay = "RNA",
    return.seurat = TRUE,
    group.by = pargs $ psbulk,
    verbose = FALSE
  )
} else {
  shared[[".psbulk"]] <- shared[["seurat"]]
}

# parse the conditions

original_ident <- Idents(shared[[".psbulk"]])

if (length(pargs $ cmp.cond) > 0 &&
      length(pargs $ with.cond) > 0) {

  cmp_mask <- parse_subset(pargs $ cmp.cond)
  with_mask <- parse_subset(pargs $ with.cond)

  # generate a temp cluster containing these files.
  intersec <- cmp_mask & with_mask
  if (sum(intersec) > 0) {
    cat(red("there should not be intersections between cmp and with."), crlf)
    stop()
  }

  temp <- rep("Others", length(original_ident))
  temp[cmp_mask] <- "Compare"
  temp[with_mask] <- "With"
  shared[[".psbulk"]][[".temp.cluster"]] <- factor(temp)
  Idents(shared[[".psbulk"]]) <- ".temp.cluster"
}

switch(
  "seurat",
  seurat = {

    if (length(pargs $ cmp.cond) > 0 &&
          length(pargs $ with.cond) > 0) {

      de <- pargs |> find_marker(
        shared[[".psbulk"]],
        ident.1 = "Compare", ident.2 = "With", verbose = FALSE
      )

      groups <- total <- shared[[".psbulk"]] |> Idents() |> levels()
      delist <- list()
      de <- de |>
        subset(p_val_adj <= pargs $ p & pct.1 > 0.10) |>
        top_n(n = pargs $ top, wt = avg_log2FC)
      genes <- rownames(de)
      colnames(de) <- c("p", "logfc", "pct1", "pct2", "p.adj")
      de $ name <- genes

      print(de |> tibble())

    } else {

      # here, we should compare each clusters separately to generate for each
      # a list of genes that are differentially expressed.

      groups <- total <- shared[[".psbulk"]] |> Idents() |> levels()
      if (pargs $ cmp.cluster[1] != "*")
        groups <- pargs $ cmp.cluster

      withclust <- pargs $ with.cluster
      if (pargs $ with.cluster[1] == "*") withclust <- NULL

      delist <- list()
      is_first <- TRUE
      for (group in groups) {
        if (pargs $ with.cluster == group) next
        if (!(group %in% names(delist))) delist[[group]] <- list()
        delist[[group]][[pargs $ with.cluster]] <- pargs |> find_marker(
          shared[[".psbulk"]],
          ident.1 = group, ident.2 = withclust, verbose = FALSE
        )

        de <- delist[[group]][[pargs $ with.cluster]] |>
          subset(p_val_adj <= pargs $ p & pct.1 > 0.10) |>
          top_n(n = pargs $ top, wt = avg_log2FC)
        genes <- rownames(de)
        colnames(de) <- c("p", "logfc", "pct1", "pct2", "p.adj")
        de $ name <- genes

        if (!is_first) cat(crlf)
        is_first = FALSE
        cat(blue(paste("differential expression between", group, 
                       "and", pargs $ with.cluster)), crlf)
        print(de |> tibble())
      }
    }

    # generating full differential expression table.

    if (stringr::str_length(pargs $ store) > 0) {
      cat(crlf)
      cat(blue("generating full differential expression table ..."), crlf)

      denames <- list()
      for (group in total) {
        for (withgroup in total) {
          if (withgroup == group) next
          if (!(group %in% names(delist))) delist[[group]] <- list()
          if (!(group %in% names(denames))) denames[[group]] <- list()
            
          cat(blue(paste("calculating differential expression between", group, 
                         "and", withgroup)), crlf)

          delist[[group]][[withgroup]] <- pargs |> find_marker(
            shared[[".psbulk"]],
            ident.1 = group, ident.2 = withgroup, verbose = FALSE
          )

          de <- delist[[group]][[withgroup]] |>
            subset(p_val_adj <= pargs $ p & pct.1 > 0.10) |>
            top_n(n = pargs $ top, wt = avg_log2FC)
          genes <- rownames(de)
          colnames(de) <- c("p", "logfc", "pct1", "pct2", "p.adj")
          de $ name <- genes

          delist[[group]][[withgroup]] <- de[, c(
            "p", "logfc", "p.adj", "name"
          )] |> tibble()
          denames[[group]][[withgroup]] <- genes
        }
      }

      shared[[pargs $ store]] <- delist
      label.markers <- lapply(denames, unlist)
      label.markers <- lapply(label.markers, unique)
      shared[[paste(pargs $ store, "names", sep = ".")]] <- label.markers
    }
  }
)

if (length(pargs $ cmp.cond) > 0 &&
      length(pargs $ with.cond) > 0) {

  Idents(shared[[".psbulk"]]) <- original_ident
  shared[[".psbulk"]][[".temp.cluster"]] <- NULL
}

shared[[".psbulk"]] <- NULL
