
parser <- argparse::ArgumentParser(
  prog = "group",
  description = "edit the grouping information sample.tsv of a run"
)

parser $ add_argument(
  "-n", dest = "novo", type = "character", default = "10x",
  help = paste(
    "indicate how to create a blueprint of grouping file when there is",
    "not. specify to 'table' to generate from the column names of a .tsv;",
    "'csv' for .csv tables, and '10x' for 10x-genomic barcodes.tsv.gz."
  )
)

parser $ add_argument(
  "-f", dest = "file", type = "character", default = "",
  help = "specify the csv or tsv matrix if you specity novo other than '10x'"
)

parser $ add_argument(
  "-r", dest = "rename", type = "character", nargs = "*", default = c(),
  help = "rename the following columns in the pre-existing sample.tsv"
)

parser $ add_argument(
  "-t", dest = "to", type = "character", nargs = "*", default = c(),
  help = "to which values should columns be renamed"
)

parser $ add_argument(
  "-c", dest = "columns", type = "character", nargs = "*", default = c(),
  help = "append columns with uniform values"
)

parser $ add_argument(
  "-v", dest = "values", type = "character", nargs = "*", default = c(),
  help = "which uniform values should the appending column take"
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

if (file.exists("sample.tsv")) {
  grp <- read.delim("sample.tsv", header = TRUE, sep = "\t",
                    comment.char = "#") |> tibble()
} else {

  if (pargs $ novo == "10x") {
    barcode <- read.delim("barcodes.tsv.gz", header = FALSE)
    grp <- data.frame(id = barcode) |> tibble()
  } else if (pargs $ novo == "table") {
    barcode <- read.delim(pargs $ file, comment.char = "#")[, -1]
    grp <- data.frame(id = colnames(barcode)) |> tibble()
  } else if (pargs $ novo == "csv") {
    barcode <- read.delim(pargs $ file, sep = ",", comment.char = "#")[, -1]
    grp <- data.frame(id = colnames(barcode)) |> tibble()
  } else {
    cat("the de-novo type:", pargs $ novo, "is not recognized.", crlf)
    stop()
  }

  colnames(grp) <- c("id")
}

if (length(pargs $ rename) != length(pargs $ to)) {
  cat("inconsistant length between -r and -t", crlf)
  stop()
}

if (length(pargs $ rename) > 0) {
  id <- 0
  for (.from in pargs $ rename) {
    id <- id + 1
    .to <- pargs $ to [id]
    if (.from %in% colnames(grp)) {
      names(grp)[names(grp) == .from] <- .to
    }
  }
}

if (length(pargs $ columns) != length(pargs $ values)) {
  cat("inconsistant length between -c and -v", crlf)
  stop()
}

if (length(pargs $ columns) > 0) {
  id <- 0
  for (.col in pargs $ columns) {
    id <- id + 1
    .val <- pargs $ values [id]
    grp[[.col]] <- rep(.val, nrow(grp))
  }
}

write.table(grp, file = "sample.tsv", quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
