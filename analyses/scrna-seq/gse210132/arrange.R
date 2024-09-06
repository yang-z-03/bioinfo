
require(fs)
require(dplyr)

series_id <- "GSE210132"

arrange_files_to_10x <- function(root_dir) {
  files <- list.files(root_dir)
  matrix <- strsplit(files, split = "_")

  i <- 1
  for (x in matrix) {
    attempt_dir <- paste(root_dir, x[1], sep = "/")
    if (!dir.exists(attempt_dir)) {
      dir.create(attempt_dir)
    }

    print(i)
    print(files[i])
    print(x[1])
    print(last(x))

    file_move(
      paste(root_dir, files[i], sep = "/"),
      paste(root_dir, x[1], last(x), sep = "/")
    )

    i <- i + 1
  }
}

root_dir <- paste("geo", series_id, sep = "/")
arrange_files_to_10x(root_dir)

sample_id_s <- function(id) {
  paste("geo/GSE210132/GSM6422", id, sep = "")
}

f <- function(id) {
  paste("geo/GSE210132/", id, sep = "/")
}

indices <- data.frame(id = 133:177)

indices $ tissue <- c(
  "circulating zsg", "tissue", "tissue", "primary", "primary", "primary",
  "metastasis", "circulating zsg", "tissue", "circulating zsg", "tissue",
  "circulating zsg", "tissue", "circulating zsg", "circulating zsg",
  "tissue", "tissue", "tissue", "circulating zsg", "circulating zsg", "tissue",
  "circulating zsg", "tissue", "circulating zsg", "tissue", "primary",
  "metastasis", "metastasis", "primary", "metastasis", "metastasis", "primary",
  "circulating zsg", "metastasis", "metastasis", "metastasis",
  "circulating zsg", "primary", "pbmc", "circulating zsg", "metastasis",
  "metastasis", "primary", "metastasis", "primary"
)

indices $ gene <- c(
  rep("wt", 5), rep("het", 8), rep("null", 3), "null", rep("null", 3),
  rep("wt", 3), "het", rep("wt", 3), rep("het", 4), "wt", "het", "wt", "het",
  "wt", "wt", "het", "wt", "het", "het", "het", "het", "het", "wt"
)

indices $ driver <- c(
  rep("prom1 > null", 3),
  "villin", "prom1", "villin", "villin",
  rep("prom1 > null", 18),
  rep("villin", 12), "prom1", "prom1 > null", rep("villin", 6)
)

indices $ tumor <- c(
  rep("no", 3),
  rep("yes", 4),
  rep("no", 18),
  rep("yes", 12), "yes", "no", rep("yes", 6)
)

indices $ loc <- c(
  "-", "-", "-", "small intestine", "gastric", "small intestine",
  "diaphragm", "-", "-", "-", "-",
  "-", "-", "-", "-",
  "-", "-", "-", "-", "-", "-",
  "-", "-", "-", "-", "small intestine",
  "mesentery", "muscle", "small intestine",
  "pancreas", "peritonium", "small intestine",
  "tumor czc", "mesentery", "ln", "ln",
  "tumor czc", "gastric", "pbmc", "tumor czc", "mesentery",
  "mesentery", "small intestine", "muscle", "small intestine"
)

saveRDS(indices, f("grouping.rds"))
readRDS(f("grouping.rds"))
