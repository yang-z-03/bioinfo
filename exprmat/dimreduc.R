
parser <- argparse::ArgumentParser(
  prog = "dimreduc",
  description = "calculate and plot dimension reduction"
)

parser $ add_argument(
  "-m", dest = "method", type = "character", default = "",
  help = paste("method of dimension reduction, accepts 'pca', 'tsne' or 'umap'")
)

parser $ add_argument(
  "-d", dest = "dimreduc", type = "character", default = "",
  help = paste("the name of existing dimension reduction to plot. specifying",
               "this tag will avoid any calculation of new reduction and only",
               "plot an existing one")
)

parser $ add_argument(
  "-r", dest = "redname", type = "character", default = "pca",
  help = paste("to which dimension reduction accepts as input")
)

parser $ add_argument(
  "--width", dest = "width", type = "double", default = 5,
  help = paste("width of the output map, in inches")
)

parser $ add_argument(
  "--height", dest = "height", type = "double", default = 5,
  help = paste("height of the output map, in inches")
)

parser $ add_argument(
  "--dim", dest = "dim", type = "integer", default = 30,
  help = paste("pca dimension components taken into calculation")
)

parser $ add_argument(
  "-t", dest = "title", type = "character", default = "",
  help = paste("the title of graphics")
)

parser $ add_argument(
  "-x", dest = "x", type = "character", default = "1",
  help = paste("x axis label")
)

parser $ add_argument(
  "-y", dest = "y", type = "character", default = "2",
  help = paste("y axis label")
)

parser $ add_argument(
  "--legend-t", dest = "ltitle", type = "character", default = "",
  help = paste("the legend (method of coloring) title")
)

parser $ add_argument(
  "--legend-l", dest = "ltext", type = "character", default = c(), nargs = "*",
  help = paste("the legend level texts, in the original order")
)

parser $ add_argument(
  "--legend-c", dest = "lcolor", type = "character", default = c(), nargs = "*",
  help = paste("the legend colors, will overwrite all groupings")
)

parser $ add_argument(
  "--draw", dest = "cond.draw", type = "character", default = c(), nargs = "*",
  help = paste("the condition sets of drawing cells")
)

parser $ add_argument(
  "--sz", dest = "size",
  type = "double", default = 0.5,
  help = paste("the point sizes of cells")
)

parser $ add_argument(
  "--highlight", dest = "cond.highlight",
  type = "character", default = c(), nargs = "*",
  help = paste("the condition sets of highlighting cells")
)

parser $ add_argument(
  "--highlight-c", dest = "highlight.color",
  type = "character", default = c("red"), nargs = "*",
  help = paste("the color sets of highlighting cells")
)

parser $ add_argument(
  "--highlight-sz", dest = "highlight.size",
  type = "double", default = c(1), nargs = "*",
  help = paste("the point sizes of highlighting cells")
)

parser $ add_argument(
  "--group-by", dest = "group.by",
  type = "character", default = "",
  help = paste("metadata column name to group by")
)

parser $ add_argument(
  "--shape-by", dest = "shape.by",
  type = "character", default = "",
  help = paste("metadata column name to identify point shapes")
)

parser $ add_argument(
  "--unmask-top", dest = "unmask.top",
  default = FALSE, action = "store_true",
  help = paste("plot the top cells last, making them stand out among others")
)

parser $ add_argument(
  "--shuffle", dest = "shuffle",
  default = FALSE, action = "store_true",
  help = paste("shuffle the locations for points too crowded")
)

parser $ add_argument(
  "--label", dest = "label",
  default = FALSE, action = "store_true",
  help = paste("label the clusters on the graph")
)

parser $ add_argument(
  "--label-sz", dest = "label.size", type = "double", default = 4,
  help = paste("the label size (if labeling on the graph)")
)

parser $ add_argument(
  "--label-c", dest = "label.color", type = "character", default = "black",
  help = paste("the label color (if labeling on the graph)")
)

parser $ add_argument(
  "--repel", dest = "repel",
  default = FALSE, action = "store_true",
  help = paste("repel the labels")
)

parser $ add_argument(
  "--alpha", dest = "alpha", type = "double", default = 1,
  help = paste("the translucency of the points")
)

parser $ add_argument(
  "--no-legend", dest = "nleg", default = FALSE,
  action = "store_true",
  help = paste("do not plot legend")
)

parser $ add_argument(
  "--rerun", dest = "rerun", default = FALSE,
  action = "store_true",
  help = paste("rerun the reduction even if there is one")
)

parser $ add_argument(
  "--tsne-perplex", dest = "perplexity", type = "integer", default = 30,
  help = paste("perplexity of t-SNE algorithm")
)

parser $ add_argument(
  "--umap-mindist", dest = "min.dist", type = "double", default = 0.3,
  help = paste("minimum distance of umap plot")
)

parser $ add_argument(
  "--umap-nn", dest = "nn", type = "integer", default = 30,
  help = paste("number of neighbors in umap manifold")
)

parser $ add_argument(
  "fname", type = "character",
  help = paste("output graphics file name")
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

if (pargs $ method == "") pargs $ method <- NULL
if (pargs $ shape.by == "") pargs $ shape.by <- NULL
if (pargs $ group.by == "") pargs $ group.by <- NULL
pargs $ ltitle <- stringr::str_replace_all(pargs $ ltitle, "_", " ")
if (pargs $ ltitle == "") pargs $ ltitle <- NULL
pargs $ title <- stringr::str_replace_all(pargs $ title, "_", " ")
if (pargs $ title == "") pargs $ title <- NULL
pargs $ x <- stringr::str_replace_all(pargs $ x, "_", " ")
pargs $ y <- stringr::str_replace_all(pargs $ y, "_", " ")

if (!shared[["is_ready"]]) {
  cat("you should normalize to generate a seurat object for later analysis")
  stop()
}

if (VariableFeatures(shared[["seurat"]]) |> is.null()) {
  shared[["seurat"]] <- FindVariableFeatures(shared[["seurat"]])
}

should_rerun <- TRUE
if (pargs $ dimreduc != "") {
  should_rerun <- FALSE
  pargs $ method <- pargs $ dimreduc # a trick, may introduce bugs
} else if (pargs $ method %in% SeuratObject::Reductions(shared[["seurat"]]) &&
             !pargs $ rerun) {
  should_rerun <- FALSE
}

if (should_rerun) {
  switch(
    pargs $ method,
    umap = {
      shared[["seurat"]] <- RunUMAP(
        shared[["seurat"]],
        reduction = pargs $ redname,
        dims = 1 : pargs $ dim,
        reduction.name = "umap",
        min.dist = pargs $ min.dist,
        n.neighbors = pargs $ nn
      )
    },
    tsne = {
      shared[["seurat"]] <- RunTSNE(
        shared[["seurat"]],
        reduction = pargs $ redname,
        dims = 1 : pargs $ dim,
        reduction.name = "tsne",
        perplexity = pargs $ perplexity
      )
    },
    pca = {
      shared[["seurat"]] <- RunPCA(
        shared[["seurat"]],
        npcs = pargs $ dim,
        reduction.name = "pca", verbose = FALSE
      )
    }
  )
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

parse_group <- function(arg) {
  l <- list()
  group <- 0
  for (item in arg) {
    if (stringr::str_starts(item, ":")) {
      group <- group + 1
      l[[group]] <- NA
    } else if (l[[group]][1] |> is.na()) {
      l[[group]] <- item
    } else {
      l[[group]] <- c(l[[group]], item)
    }
  }

  r <- list()
  group <- 0
  cellnames <- Cells(shared[["seurat"]])
  for (item in l) {
    group <- group + 1
    r[[group]] <- cellnames[parse_subset(l[[group]])]
  }

  return(r)
}

drawnames <- NULL
if (pargs $ cond.draw |> length() > 0)
  drawnames <- parse_subset(pargs $ cond.draw)

highlights <- NULL
ontop <- NULL
if (pargs $ cond.highlight |> length() > 0) {
  highlights <- parse_group(pargs $ cond.highlight)
  if (pargs $ unmask.top) ontop <- purrr::reduce(highlights, c)
}

# this method destruct the original seurat object.

#. if (pargs $ ltext |> length() > 0) {
#.   new_ident <- setNames(pargs $ ltext, levels(shared[["seurat"]]))
#.   shared[["seurat"]] <- RenameIdents(shared[["seurat"]], new_ident)
#. }

if (highlights |> is.null()) {
  dimplot <- Seurat::DimPlot(
    object = shared[["seurat"]],
    cells = drawnames,
    pt.size = pargs $ size,
    reduction = pargs $ method,
    group.by = pargs $ group.by,
    shape.by = pargs $ shape.by,
    shuffle = pargs $ shuffle,
    label = pargs $ label,
    label.size = pargs $ label.size,
    label.color = pargs $ label.color,
    repel = pargs $ repel,
    alpha = pargs $ alpha
  )
} else {
  dimplot <- Seurat::DimPlot(
    object = shared[["seurat"]],
    cells = drawnames,
    pt.size = pargs $ size,
    reduction = pargs $ method,
    group.by = pargs $ group.by,
    shape.by = pargs $ shape.by,
    shuffle = pargs $ shuffle,
    label = pargs $ label,
    label.size = pargs $ label.size,
    label.color = pargs $ label.color,
    repel = pargs $ repel,
    alpha = pargs $ alpha,
    cells.highlight = highlights,
    cols.highlight = pargs $ highlight.color,
    sizes.highlight = pargs $ highlight.size,
    order = ontop
  )
}

if (!(pargs $ title |> is.null())) {
  dimplot <- dimplot + labs(title = pargs $ title)
}

if (!(pargs $ ltitle |> is.null())) {
  dimplot <- dimplot + labs(color = pargs $ ltitle)
}

if (pargs $ nleg) {
  dimplot <- dimplot + NoLegend()
}

if (pargs $ ltext |> length() > 0) {
  dimplot <- dimplot +
    ggplot2::scale_color_manual(values = pargs $ lcolor, labels = pargs $ ltext)
}

dimplot <- dimplot +
  labs(x = pargs $ x, y = pargs $ y) +
  unify_theme_font()

ggsave(
  pargs $ fname, dimplot, units = "in",
  width = pargs $ width, height = pargs $ height, dpi = 600
)
