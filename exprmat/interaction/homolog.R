
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

key <- convtable[convtable $ entrez1 %in% entrezids, ]
key <- key[key $ tax2 == 9606, ]
key2 <- convtable[convtable $ entrez2 %in% entrezids, ]
key2 <- key2[key2 $ tax1 == 9606, ]
key2 <- data.frame(
  tax1 = key2 $ tax2,
  entrez1 = key2 $ entrez2,
  relation = key2 $ relation,
  tax2 = key2 $ tax1,
  entrez2 = key2 $ entrez1
)

# remove duplicates within key1 or key2: we should just abandon all of the
# copies. since multi-mapping orthologs are directly ignored.

has_dupl <- key $ entrez2[key $ entrez2 |> duplicated()] |> unique()
key <- key[!(key $ entrez2 %in% has_dupl), ]
has_dupl <- key2 $ entrez2[key2 $ entrez2 |> duplicated()] |> unique()
key2 <- key2[!(key2 $ entrez2 %in% has_dupl), ]

key <- rbind(key, key2) |> tibble()
key <- key[!duplicated(key), ]

cat(blue("ortholog correspondence table:", crlf, crlf))
print(key)

orthoinfo <- match(key $ entrez2, targetdb $ entrez)
flag <- !is.na(orthoinfo)
orthoinfo <- orthoinfo[flag]
key <- key[flag, ]

cat(crlf)
cat(blue("one -> one ortholog table"), crlf)
print(targetdb[orthoinfo, ] |> tibble())

targetortho <- targetdb[orthoinfo, ]
matches <- match(shared $ meta_gene $ entrez, key $ entrez1)
orthologs <- targetortho[matches, ] $ symbol
human_entrez <- targetortho[matches, ] $ entrez

orthotable <- data.frame(
  gene = shared $ meta_gene $ gene,
  entrez = shared $ meta_gene $ entrez,
  ortholog = orthologs,
  entrez.human = human_entrez
)

# save the files

if (!dir.exists("cpdb")) dir.create("cpdb")

cat(crlf)
cat(blue("writing homolog table ..."), crlf)
write.table(orthotable, file = "cpdb/orthologs.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
