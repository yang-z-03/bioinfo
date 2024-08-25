
library(magrittr)

# The following list indicates the code and description of the issues detected
# in GFF3 files:
#
# - `NCOLUMNS_EXCEEDED`: Input file contains lines with more than 9 fields.
# - `NCOLUMNS_INFERIOR`: Input file contains lines with less than 9 fields
# - `TOO_MANY_FEATURE_TYPES`: Input file contains too many (more than 100)
#                             different feature types.
# - `NO_ID`: ID attribute not found in any feature
# - `DUPLICATED_ID`: There are duplicated IDs
# - `ID_IN_MULTIPLE_CHR`: The same ID has been found in more than one chromosome
# - `NO_PARENT`: Parent attribute not found in any feature
# - `MISSING_PARENT_ID`: There are missing Parent IDs
# - `PARENT_INDIFFERENTCHR`: There are features whose parent is located in a
#                            different chromosome.
# - `PARENT_DEFINED_BEFORE_ID`: Feature ids referenced in Parent attribute
#                               before being defined as ID.
# - `NOT_GROUPED_BY_CHR`: Features are not grouped by chromosome
# - `NOT_SORTED_BY_COORDINATE`: Features are not sorted by start coordinate
# - `NOT_VALID_WARNING`: File cannot be recognized as valid GFF3.
# - `NOT_VALID_ERROR`: File cannot be recognized as valid GFF3. Parsing errors.
#
# The following list indicates the code and description of the issues detected
# in GTF files:
#
# - `NCOLUMNS_EXCEEDED`: Input file contains lines with more than 9 fields
# - `NCOLUMNS_INFERIOR`: Input file contains lines with less than 9 fields
# - `TOO_MANY_FEATURE_TYPES`: Input file contains too many (more than 100)
#                             different feature types
# - `NO_GENE_ID_ATTRIBUTE`: gene_id attribute not found in any feature
# - `MISSING_GENE_IDs`: There are features without gene_id attribute
# - `NO_GENE_FEATURES`: Gene features are not included in this GTF file
# - `DUPLICATED_GENE_IDs`: There are duplicated gene_ids
# - `GENE_ID_IN_MULTIPLE_CHR`: The same gene_id has been found in more than one
#                              chromosome
# - `NO_TRANSCRIPT_ID_ATTRIBUTE`: transcript_id attribute not found in any
#                                 feature There are no elements with
#                                 `transcript_id` attribute
# - `MISSING_TRANSCRIPT_IDs`: There are features without transcript_id attribute
# - `NO_TRANSCRIPT_FEATURES`: Transcript features are not included in this file
# - `DUPLICATED_TRANSCRIPT_IDs`: There are  duplicated transcript_ids
# - `TRANSCRIPT_ID_IN_MULTIPLE_CHR`: The same transcript_id has been found in
#                                    more than one chromosome
# - `DUPLICATED_GENE_AND_TRANSCRIPT_IDs`: Same id has been defined as `gene_id`
#                                         and `transcript_id`
# - `NOT_GROUPED_BY_CHR`: Features are not grouped by chromosome
# - `NOT_SORTED_BY_COORDINATE`: Features are not sorted by start coordinate
# - `NOT_VALID_WARNING`: File cannot be recognized as valid GTF.
# - `NOT_VALID_ERROR`: File cannot be recognized as valid GTF. Parsing errors.
#
# returns:  A data frame of detected issues, including a short code name,
#           a description and estimated severity each. If no issues are detected
#           the function will return an empty data frame.
#
gff3_check <- function(gff_file, gff_table = gff3_list(gff_file)) { # nolint

  if (!base::file.exists(gff_file)) {
    stop(paste0("Input GFF3 file ", gff_file, " not found."))
  }

  options(dplyr.summarise.inform = FALSE)

  errors_det <- data.frame(error = character(),
                           message = character(),
                           severity = character())

  detect_errors <- function() {
    fields_n <- table(utils::count.fields(gff_file, sep = "\t",
                                          comment.char = "#", quote = ""))

    if (max(as.numeric(names(fields_n))) > 9) {
      # exceeds 9 columns
      err_msg <- paste0("Input file contains lines with more than 9 fields")
      message(err_msg)
      errors_det[nrow(errors_det) + 1, ] <-
        c("NCOLUMNS_EXCEEDED", err_msg, "HIGH")

    } else if (min(as.numeric(names(fields_n))) < 9) {
      # less than 9 columns
      err_msg <- paste0("Input file contains lines with less than 9 fields")
      message(err_msg)
      errors_det[nrow(errors_det) + 1, ] <-
        c("NCOLUMNS_INFERIOR", err_msg, "HIGH")
    }

    table_data <- gff_table

    # feature types too much
    different_features <- unique(table_data $ feature)
    if (length(different_features) > 100) {
      err_msg <- paste0("There are too many different feature types (",
                        length(different_features), ")")
      message(err_msg)
      errors_det[nrow(errors_det) + 1, ] <-
        c("TOO_MANY_FEATURE_TYPES", err_msg, "MEDIUM")
    }

    id_pattern <- "^(.*;)?( *ID *= *([^;]+))(.*)$"
    selected_id_lines <- grep(id_pattern, table_data $ attribute,
                              perl = TRUE, ignore.case = TRUE)
    if (length(selected_id_lines) == 0) {
      err_msg <- paste0("ID attribute not found in any feature")
      message(err_msg)
      errors_det[nrow(errors_det) + 1, ] <- c("NO_IDs", err_msg, "HIGH")

    } else {
      id_df <-
        data.frame(ChrID = table_data[selected_id_lines, "seqname"],
                   FeatureID = table_data[selected_id_lines, "feature"],
                   ID = gsub(id_pattern, "\\3",
                             table_data[selected_id_lines, "attribute"],
                             perl = TRUE, ignore.case = TRUE),
                   PosID = selected_id_lines)

      # check repeated ids for different feature types.

      id_cnt <- id_df %>%
        dplyr::group_by(.data $ FeatureID, .data $ ID) %>%
        dplyr::summarise(N = dplyr::n(), FirstID = min(.data$PosID))

      if (length(id_cnt$ID) != length(unique(id_cnt$ID))) {
        err_msg <- paste0("There are ",
                          length(id_cnt$ID) - length(unique(id_cnt$ID)),
                          " duplicated IDs")
        message(err_msg)
        errors_det[nrow(errors_det) + 1, ] <-
          c("DUPLICATED_IDs", err_msg, "HIGH")
      }

      if (any(duplicated(unique(id_df[, c("ID", "ChrID")])$ID))) {
        err_msg <- paste0("The same  ", length(which(duplicated(unique(
          id_df[, c("ID", "ChrID")])$ID))),
          " ids have been found in more than one chromosome")
        message(err_msg)
        errors_det[nrow(errors_det) + 1, ] <-
          c("ID_IN_MULTIPLE_CHR", err_msg, "HIGH")
      }
    }

    parent_pattern <- "^(.*;)?( *Parent *= *([^;]+))(.*)$"
    selected_parent_lines <- grep(parent_pattern, table_data $ attribute,
                                  perl = TRUE, ignore.case = TRUE)

    if (length(selected_parent_lines) == 0) {
      err_msg <- paste0("Parent attribute not found in any feature")
      message(err_msg)
      errors_det[nrow(errors_det) + 1, ] <- c("NO_PARENTs", err_msg, "MEDIUM")

    } else {
      parent_df <-
        data.frame(ChrChild = table_data[selected_parent_lines, "seqname"],
                   Feature = table_data[selected_parent_lines, "feature"],
                   ID = gsub(parent_pattern, "\\3",
                             table_data[selected_parent_lines, "attribute"],
                             perl = TRUE, ignore.case = TRUE),
                   PosChild = selected_parent_lines)

      # multiple parent occurances.
      multi_parent <- grep(",", parent_df $ ID)
      if (length(multi_parent) > 0) {
        multiparent_df <-
          as.data.frame(do.call("rbind", lapply(multi_parent, function(x) {
            t(sapply(strsplit(parent_df $ ID[x], ",")[[1]], function(y) {
              c(parent_df $ ChrChild[x],
                parent_df $ Feature[x], y,
                as.numeric(parent_df $ PosChild[x]))
            }))
          })))
        colnames(multiparent_df) <- c("ChrChild", "Feature", "ID", "PosChild")
        multiparent_df $ PosChild <- as.numeric(multiparent_df $ PosChild)
        parent_df <- rbind(parent_df[-multi_parent, ], multiparent_df)
      }

      parent_df <- parent_df %>%
        dplyr::group_by(.data $ Feature, .data $ ID, .data $ ChrChild) %>%
        dplyr::summarise(n_children = dplyr::n(),
                         first_child = min(.data $ PosChild)) %>%
        dplyr::arrange(.data $ ID)

      # with missing parents
      merged_df <- merge(x = parent_df, y = id_df, by = "ID", all.x = TRUE)
      missing_df <- which(is.na(merged_df $ FeatureID))
      if (length(missing_df) > 0) {
        err_msg <- paste0("There are ", length(missing_df),
                          " missing Parent IDs")
        message(err_msg)
        errors_det[nrow(errors_det) + 1, ] <-
          c("MISSING_PARENT_IDs", err_msg, "HIGH")
      }

      if (any(merged_df $ ChrChild != merged_df $ ChrID, na.rm = TRUE)) {
        err_msg <-
          paste0("There are ",
                 length(which(merged_df$ChrChild != merged_df$ChrID)),
                 " features whose Parent located in a different chromosome")
        message(err_msg)
        errors_det[nrow(errors_det) + 1, ] <-
          c("PARENT_IN DIFFERENT CHR", err_msg, "HIGH")
      }

      # feature ids referenced in Parent attribute before being defined
      refer_before_defined <-
        which(merged_df $ first_child - merged_df $ Pos < 0)
      if (length(refer_before_defined) > 0) {
        err_msg <- paste0(length(refer_before_defined),
                          " feature ids referenced in Parent attribute" +
                            " before being defined as ID")
        message(err_msg)
        errors_det[nrow(errors_det) + 1, ] <-
          c("PARENT_DEFINED_BEFORE_ID", err_msg, "MEDIUM")
      }
    }

    # chromosome order
    chr_order <- unique(table_data $ seqname)
    chr_factor <- factor(table_data $ seqname, levels = chr_order)
    if (is.unsorted(chr_factor)) {
      err_msg <- paste0("Features are not grouped by chromosome")
      message(err_msg)
      errors_det[nrow(errors_det) + 1, ] <-
        c("NOT_GROUPED_BY_CHR", err_msg, "MEDIUM")
    }

    # start coordination order
    is_start_sorted <- table_data %>%
      dplyr::group_by(.data $ seqname) %>%
      dplyr::summarise(unsort = is.unsorted(.data $ start, strictly = FALSE))

    if (any(is_start_sorted $ unsort)) {
      err_msg <- paste0("Features are not sorted by start coordinate in ",
                        sum(is_start_sorted $ unsort), " chromosomes")
      message(err_msg)
      errors_det[nrow(errors_det) + 1, ] <-
        c("NOT_SORTED_BY_COORDINATE", err_msg, "MEDIUM")
    }

    if (nrow(errors_det) == 0) {
      message(paste0(gff_file, " OK: No errors detected."))
    }
    return(errors_det)
  }

  errors_det <- detect_errors()

  rm(list = setdiff(ls(), "errors_det"))
  gc(reset = TRUE)
  return(errors_det)
}
