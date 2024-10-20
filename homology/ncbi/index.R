
library(tibble)
library(dplyr)

info <- read.table(
  "homology/ncbi/ncbi-genes.gz", header = FALSE,
  sep = "\t", comment.char = "#", quote = "\"", fill = TRUE
)

colnames(info) <- c(
  "taxid", "entrez", "symbol", "locus", "synonyms", "dbxref", "chr", "map",
  "description", "type", "symbol.auth", "name.auth", "nomenclature", "other",
  "date", "feature"
)

taxtable <- info $ taxid |> table()
dir.create("homology/ncbi/genes")
keys <- names(taxtable)
keys <- keys[taxtable > 1000]

# we only export organisms with greater than 1000 genes.
# since we currently do not care about those rare species and prokaryotes.

n <- 1
for (i in keys) {
  species <- info[info $ taxid == as.integer(i), ]
  cat("saving species", i, "(total of", nrow(species), "genes)",
      "-", n, "/", length(taxtable), "\n")
  saveRDS(species, paste("homology/ncbi/genes/", i, ".rds", sep = ""))
  n <- n + 1
}
