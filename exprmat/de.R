
parser <- argparse::ArgumentParser(
  prog = "de",
  description = "get a list of differential expression genes amoung clusters"
)

parser $ add_argument(
  "-m", dest = "method", type = "character", default = "",
  help = paste("method of dimension reduction, accepts 'singler' or 'seurat'")
)

parser $ add_argument(
  "-d", dest = "direc", type = "character", default = "",
  help = paste("the direction of expression change, one of 'up' and 'down'")
)

parser $ add_argument(
  "--cmp", dest = "cmp.cluster",
  type = "character", default = c(), nargs = "+",
  help = paste("cluster names (combined with OR) as the comparing set")
)

parser $ add_argument(
  "--with", dest = "with.cluster",
  type = "character", default = c(), nargs = "+",
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

original_ident <- shared[["seurat"]] @ active.ident
if (length(pargs $ cmp.cond) > 0 &&
      length(pargs $ with.cond) > 0) {

}