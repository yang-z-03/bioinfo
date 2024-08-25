
library(magrittr)
library(tibble)

gff3_summary <- function(gff_file) {

  n_fields <- max(utils::count.fields(gff_file, sep = "\t",
                                      comment.char = "#", quote = ""))

  if (n_fields > 9) {
    message("[warning] input file contains lines with more than 9 fields")
  }

  gff_table <-
    utils::read.table(gff_file, sep = "\t", comment.char = "#",
                      colClasses = c("character", "NULL", "character",
                                     "numeric", "numeric", "NULL",
                                     "NULL", "NULL", "NULL"),
                      stringsAsFactors = FALSE, quote = "")
  colnames(gff_table) <- c("seqname", "feature", "start", "end")

  return(gff3_summary_from_table(gff_table))
}

gff3_summary_from_table <- function(gff_table) {
  gff_summ <- gff_table %>%
    dplyr::group_by(.data $ feature) %>%
    dplyr::summarise(average_length = mean(.data $ end - .data $ start),
                     max_length = max(.data $ end - .data $ start),
                     min_length = min(.data $ end - .data $ start),
                     n = dplyr::n())

  if (!("chromosome" %in% gff_summ $ feature)) {
    gff_summ <- gff_summ %>%
      tibble::add_row(type = "chromosome",
                      n = length(unique(gff_table $ seqname)), .before = 1)
  }

  rm(list = setdiff(ls(), "gff_summ"))
  gc(reset = TRUE)
  return(gff_summ)
}