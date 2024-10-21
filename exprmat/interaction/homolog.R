
parser <- argparse::ArgumentParser(
  prog = "cpdbhom",
  description = "dump a homolog table for cpdb"
)

parser $ add_argument(
  "-s", dest = "species", type = "integer", default = 9606,
  help = paste("which species does this dataset comes from")
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

# the program will automatically dumps a table contains the genes and 
# corresponding homologs from the project. it will by default converts
# to human homologs, because cellphonedb is a human database.

if (pargs $ species == 9606) {
  cat('the source species is already human!', crlf)
  cat('no operation performed', crlf)
  stop()
}

# query the homologs

if (! file.exists("~/bioinfo/homology/ncbi/orthologs.rds")) {
  convtable <- read.table(
    "~/bioinfo/homology/ncbi/orthologs.gz", header = FALSE,
    sep = "\t", comment.char = "#", quote = "\"", fill = TRUE
  )

  colnames(convtable) <- c(
    "tax1", "entrez1", "relation", "tax2", "entrez2"
  )
  saveRDS(convtable, "~/bioinfo/homology/ncbi/orthologs.rds")
} else convtable <- readRDS("~/bioinfo/homology/ncbi/orthologs.rds")

entrezids <- shared $ meta_gene $ entrez
targetdb <- paste("~/bioinfo/homology/ncbi/genes/9606.rds", sep = "")
targetdb <- readRDS(targetdb)

# the entrezids to be queried may rest in either entrez1 or entrez2.
# look for convtable to find them.

orthologs <- NULL
human_entrez <- NULL

for (id in entrezids) {
  key <- convtable[convtable $ entrez1 == id, ]
  key <- key[key $ tax2 == pargs $ target, ]
  key2 <- convtable[convtable $ entrez2 == id, ]
  key2 <- key2[key2 $ tax1 == pargs $ target, ]
  key2 <- data.frame(
    tax1 = key2 $ tax2,
    entrez1 = key2 $ entrez2,
    relation = key2 $ relation,
    tax2 = key2 $ tax1,
    entrez2 = key2 $ entrez1
  )

  key <- rbind(key, key2) |> tibble()
  key <- key[!duplicated(key), ]

  orthoinfo <- match(key $ entrez2, targetdb $ entrez)
  flag <- !is.na(orthoinfo)
  orthoinfo <- orthoinfo[flag]
  key <- key[flag, ]

  if (length(orthoinfo) == 1) {
    orthologs <- c(orthologs, ortholog = targetdb[orthoinfo, ] $ symbol)
    human_entrez <- c(orthologs, ortholog = targetdb[orthoinfo, ] $ entrez)
  } else if (length(orthoinfo) > 1) {
    orthologs <- c(orthologs, ".multi")
    human_entrez <- c(human_entrez, ".multi")
  } else {
    orthologs <- c(orthologs, ".none")
    human_entrez <- c(human_entrez, ".none")
  }
}

orthotable <- data.frame(
  gene = shared $ meta_gene $ gene,
  entrez = shared $ meta_gene $ entrez,
  ortholog = orthologs,
  entrez.human = human_entrez
)

# save the files

if (!dir.exists("cpdb")) dir.create("cpdb")

cat(blue("writing homolog table ..."), crlf)
write.table(orthotable, file = "cpdb/orthologs.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
