
library(data.tree)
library(data.table)
library(magrittr)
library(RJSONIO)

gff3_features <- function(gff_file, include_counts = FALSE, # nolint
                          out_format = c("tree", "data.frame", "json")) {

  out_format <- match.arg(out_format)

  pair_file <- paste0(gff_file, ".pairs")
  if (!base::file.exists(pair_file)) {
    message("creating pairs file ...")
    pair_file <- gff3_pairs(gff_file, pair_file)
  }

  if (is.data.frame(pair_file)) {
    pair_table <- pair_file
  } else if (is.vector(pair_file) &&
               is.character(pair_file) &&
               length(pair_file) == 1) {
    if (file.exists(pair_file)) {
      pair_table <-
        utils::read.table(pair_file, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE,
                          colClasses = c("character", "character", "numeric"))
    } else {
      stop("only data.frame of pairs or a valid path of pairs file are allowed")
    }
  }

  if (!all(c("element", "parent") %in% names(pair_table))) {
    stop("missing element or parent columns in pair data")
  }

  pair_table <- pair_table %>%
    dplyr::relocate(.data $ parent, .before = .data $ element)

  data_tree <- data.tree::as.Node(pair_table, mode = "network")
  if (out_format == "tree") {
    return(data_tree)

  } else {
    if (include_counts) {
      descendants <- function(node) {
        descends <-
          data.tree::ToDataFrameTree(node, "name", "n")[-1, c("name", "n")]
        aggr_data <- descends %>%
          dplyr::group_by(.data $ name) %>%
          dplyr::summarise(n = sum(.data $ n))
        return(paste(paste(aggr_data $ name,
                           aggr_data $ n, sep = ":"),
                     collapse = " "))
      }

    } else {
      descendants <- function(node) {
        descends <-
          unique(data.tree::ToDataFrameTree(node, "name")[-1, c("name")])
        return(paste(descends, collapse = " "))
      }
    }

    data_tree $ Do(function(node) node $ descendants <- descendants(node))

    q_features <- as.data.frame((data_tree$Get("descendants"))[-1])
    colnames(q_features) <- c("blocks")

    if (include_counts) {
      qnames <-
        dplyr::left_join(
          data.frame(name = names((data_tree $ Get("descendants"))[-1])), # nolint
                     data.tree::ToDataFrameTree(data_tree, "name", "n")[-1, c("name", "n")] %>% # nolint
                       dplyr::group_by(.data $ name) %>%
                       dplyr::summarise(n = sum(.data $ n)))
      q_features $ features <- with(qnames, paste0(name, ":", n))

    } else {
      q_features $ features <- names((data_tree $ Get("descendants"))[-1])
    }

    q_features <- unique(q_features)

    if (out_format == "data.frame") {
      return(q_features)
    } else if (out_format == "json") {
      feature_list <- as.list(as.data.frame(t(q_features)))
      json <- RJSONIO::toJSON(list(features = feature_list),
                              .escapeEscapes = FALSE)

      return(json)
    }
  }
}