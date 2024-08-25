
library(magrittr)
library(data.table)
library(withr)

gff3_pairs <- function(gff_file, pair_file,
                       force_overwrite = FALSE) {

  if (!base::file.exists(gff_file)) {
    stop(paste0("Input GFF3 file ", gff_file, " not found."))
  }

  if (missing(pair_file)) {
    foutput <- paste0(gff_file, ".pairs")
  } else {
    foutput <- pair_file
  }

  if (base::file.exists(foutput) && force_overwrite == FALSE) {
    if (!missing(pair_file)) {
      message("Pairs file for this GFF already exists.")
    }

    pair_table <-
      utils::read.table(foutput, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE,
                        colClasses = c("character", "character", "numeric")) |>
      as.data.table()
    return(pair_table)
  }

  gff_table <- gff3_list(gff_file)

  id_pattern <- "^(.*;)?( *id *= *([^;]+))(.*)$"
  selected_id_lines <- grep(id_pattern, gff_table $ attribute,
                            perl = TRUE, ignore.case = TRUE)

  has_id <- data.table(feature_id = gff_table[selected_id_lines, "feature"],
                       id = gsub(id_pattern, "\\3",
                                 gff_table[selected_id_lines, "attribute"],
                                 perl = TRUE, ignore.case = TRUE))

  parent_pattern <- "^(.*;)?( *Parent *= *([^;]+))(.*)$"
  selected_parent_lines <- grep(parent_pattern, gff_table$attribute,
                                perl = TRUE, ignore.case = TRUE)

  has_parent <-
    data.table(feature = gff_table[selected_parent_lines, "feature"],
               id = gsub(parent_pattern, "\\3",
                         gff_table[selected_parent_lines, "attribute"],
                         perl = TRUE, ignore.case = TRUE))

  # multiple parent occurences
  with_multiple_parents <- grep(",", has_parent $ id)
  if (length(with_multiple_parents) > 0) {
    items <- do.call("rbind", lapply(with_multiple_parents, function(x) {
      t(sapply(strsplit(has_parent $ id[x], ",")[[1]], function(y) {
        c(has_parent $ feature[x], y)
      }))
    }))
    colnames(items) <- c("feature", "id")
    has_parent <- rbind(has_parent[-with_multiple_parents, ], items)
  }

  merge_df <- merge(x = has_parent, y = has_id, by = "id", all.x = TRUE)
  count_pairs <- merge_df %>% dplyr::count(.data $ feature, .data $ feature_id)

  # missing parent ids
  if (any(is.na(count_pairs $ feature_id))) {

    message("warning: there are missing parent IDs and will be ignored.")
    # Option 1: (active) For features with missing parents, parent feature
    # changed from NA to empty string in pairs

    count_pairs[is.na(count_pairs$feature_id), ] $ feature_id <- ""
    # Option 2: Features with missing parents are removed from pairs file
    # count_pairs <- count_pairs[!is.na(count_pairs $ feature_id),] # nolint
  }

  element_counts <-
    data.table(elt = gff_table[-selected_parent_lines, "feature"],
               p = rep("", nrow(gff_table) - length(selected_parent_lines))) %>%
    dplyr::count(.data $ elt, .data $ p)

  colnames(count_pairs) <- c("element", "parent", "n")
  colnames(element_counts) <- c("element", "parent", "n")

  df <- rbind(count_pairs, element_counts)
  withr::with_options(c(scipen = 999),
                      utils::write.table(df, file = foutput, col.names = TRUE,
                                         row.names = FALSE, quote = FALSE,
                                         sep = "\t", append = FALSE))

  rm(list = setdiff(ls(), "df"))
  gc(reset = TRUE)
  return(df)
}