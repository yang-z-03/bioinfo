
library(magrittr)
library(biomartr)
library(stringr)

# refseq ----------------------------------------------------------------------

# read the organism list exist locally
local <- readRDS("refseq/organisms.rds")

# the desired list, change this list if you want to add databases and run the
# script from the project directory
wanted <- readRDS("refseq/wanted.rds")

names <- readRDS("tax/ncbi-names.rds")
scientific_names <- names[names$name_class == "scientific name", , drop = FALSE]
scientific_names <- scientific_names[scientific_names$name_txt %in% wanted, ]

rm(names)
gc()

requested <- scientific_names[!scientific_names["tax_id"] %in% local]

for (id in requested$tax_id) {

  if (id %in% local) {
    cat("The taxon genome", id, "exists locally.\n")
    next
  }

  try({

    # download the sequence from refseq
    path <- biomartr::getCollection(db = "refseq",
                                    organism = id,
                                    reference = TRUE,
                                    analyse_genome = FALSE,
                                    path = "d:/projects/r")

    # if downloaded successfully, append to the downloaded list
    # if (str_length(path) > 0) local <- c(local, id) # nolint
  })
}

saveRDS(local, file = "refseq/organisms.rds")

# ensembl ---------------------------------------------------------------------

local <- readRDS("ensembl/organisms.rds")
wanted <- readRDS("ensembl/wanted.rds")
requested <- scientific_names[!scientific_names["tax_id"] %in% local]

for (id in requested$tax_id) {

  if (id %in% local) {
    cat("The taxon genome (ensembl)", id, "exists locally.\n")
    next
  }

  try({

    # download the sequence from refseq
    path <- biomartr::getCollection(db = "ensembl",
                                    organism = id,
                                    reference = TRUE,
                                    analyse_genome = FALSE,
                                    path = "d:/projects/r")

    # if downloaded successfully, append to the downloaded list
    # if (str_length(path) > 0) local <- c(local, id) # nolint
  })
}
