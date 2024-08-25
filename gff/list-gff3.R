
library(data.table)

gff3_list <- function(gff_file) {

  n_fields <- max(utils::count.fields(gff_file, sep = "\t",
                                      comment.char = "#", quote = ""))

  if (n_fields > 9) {
    message("[warning] input file contains lines with more than 9 fields")
  }

  cnames <- c("seqname", "source", "feature", "start", "end", "score",
              "strand", "frame", "attribute")
  table_data <- utils::read.delim(gff_file, header = FALSE, sep = "\t",
                                  comment.char = "#", quote = "",
                                  blank.lines.skip = TRUE) |>
    `colnames<-`(cnames) |>
    as.data.table()

  return(table_data)
}

gff3_list_attributes <- function(gff_table, filter_feature, patterns,
                                 column_names) {

  gff_selected <- gff_table[feature == filter_feature] # nolint
  nrow_gff <- nrow(gff_selected)

  list <- lapply(patterns, function(pat) {
    pattern <- paste0("^(.*;)?( *", pat, " *= *([^;]+))(.*)$")
    selected_lines <- grep(pattern, gff_selected $ attribute,
                           perl = TRUE, ignore.case = TRUE)
    empty <- data.table()
    empty[, pat] <- rep("", nrow_gff)
    empty[, ".by"] <- 1:nrow_gff
    empty[selected_lines, pat] <-
      gsub(pattern, "\\3", gff_selected[selected_lines, ] $ attribute,
           perl = TRUE, ignore.case = TRUE)

    return(empty)
  })

  l <- Reduce(merge, list)[, -1]
  colnames(l) <- column_names
  return(cbind(gff_selected[, -9], l))
}