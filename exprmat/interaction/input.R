
parser <- argparse::ArgumentParser(
  prog = "cpdbin",
  description = "generate the inputs required by cellphonedb"
)

parser $ add_argument(
  "--subset", dest = "cond",
  type = "character", default = c(), nargs = "+",
  help = paste("conditions specifying the cluster and groupings")
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

filter <- parse_subset(l[[item]])
filter_expr <- shared[["seurat"]][, filter]
filter_cell <- shared[["meta_sample"]][, filter]

# the cell names are auto-generated.