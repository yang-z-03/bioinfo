
library(magrittr)
library(data.table)
library(withr)

gff3_sort <- function(gff_file, sorted_output_file, # nolint
                      force_overwrite = FALSE) {

  if (!base::file.exists(gff_file)) {
    stop(paste0("Input GFF3 file ", gff_file, " not found."))
  }

  if (missing(sorted_output_file)) {
    foutput <- paste0(sub(pattern = "(.*)\\..{0,4}$",
                          replacement = "\\1", gff_file), ".sorted.gff3")
  } else {
    foutput <- sorted_output_file
  }

  if (base::file.exists(foutput) && force_overwrite == FALSE) {
    message("sorted file for this gff already exists.")
    return(foutput)
  }

  file_connection <- file(gff_file, "r")
  gff_header <- c()
  read_header <- TRUE
  comment_pattern <- "(^#(.*)$)|(^#([\t ]*)$)"
  blank_pattern <- "(^[\t ]*$)"

  # find and read the header line.
  while (read_header) {
    current_line <- readLines(file_connection, n = 1)
    if (identical(current_line, character(0))) {
      read_header <- FALSE
    } else if (grepl(comment_pattern, current_line, perl = TRUE)) {
      gff_header <- c(gff_header, current_line)
    } else if (grepl(blank_pattern, current_line, perl = TRUE)) {
    } else {
      read_header <- FALSE
    }
  }

  close(file_connection)

  tree <- gff3_features(gff_file, out_format = "tree")

  feature_list <- rev(unique(rev((sapply(
    data.tree::Traverse(tree, traversal = "level"), # nolint
    function(node) {{ node$name }})[-1]))))

  # feature_order <- factor(feature_list, levels = feature_list) # nolint
  gff_table <- utils::read.delim(gff_file, comment.char = "#", header = FALSE,
                                 sep = "\t", blank.lines.skip = TRUE)

  gff_ordered <- gff_table %>%
    dplyr::arrange(.data $ V1, .data $ V4, dplyr::desc(.data$V5),
                   factor(.data$V3, levels = feature_list))

  file_connection <- file(foutput, "w+")

  if (length(gff_header) > 0) {
    writeLines(gff_header, file_connection)
  }

  withr::with_options(c(scipen = 999),
                      utils::write.table(gff_ordered, file_connection,
                                         append = TRUE, sep = "\t",
                                         quote = FALSE, row.names = FALSE,
                                         col.names = FALSE))

  close(file_connection)

  rm(list = setdiff(ls(), "foutput"))
  gc(reset = TRUE)
  return(foutput)
}
