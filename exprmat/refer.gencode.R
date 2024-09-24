
parser <- argparse::ArgumentParser(
  prog = "refer.gencode",
  description = "select an organism as reference genome using gencode library"
)

parser $ add_argument(
  type = "character", dest = "taxo",
  help = "select reference genome for given taxo"
)

parser $ add_argument(
  "-a", type = "character", dest = "annot",
  help = "the appending annotation hub index"
)

parser $ add_argument(
  "-s", "--search", action = "store_const",
  help = "search in the annotation hub for valid annotations",
  dest = "method", const = "search", default = "construct"
)

parser $ add_argument(
  "-q", dest = "query", type = "character", nargs = "*", default = c(),
  help = "query vector when searching with --search"
)

parser $ add_argument(
  "-r", "--rename", action = "store_true",
  help = "ask user to rename the column before saving",
  dest = "rename", default = FALSE
)

parser $ add_argument(
  "-m", dest = "mt", type = "character", default = "chrM",
  help = "the sequence name of the mitochondrial genome"
)

parser $ add_argument(
  "--s-phase", dest = "sgenes", type = "character", nargs = "*", default = c(),
  help = "specify the s phase markers for current taxo"
)

parser $ add_argument(
  "--g2m-phase", dest = "g2mgenes", type = "character",
  nargs = "*", default = c(),
  help = "specify the g-to-m phase markers for current taxo"
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

# search a given taxo name in the annotation hub.
# these annotations are largely based on ensembl. while we choose entrez as
# the primary annotation keys since they are more beautiful :)

if (pargs $ method == "search") {

  suppressPackageStartupMessages(
    require(AnnotationHub)
  )

  hub <- AnnotationHub()
  result <- query(hub, c("orgdb", pargs $ query))
  rawtable <- mcols(result)
  rnames <- rownames(rawtable)
  table <- rawtable |> as.data.frame() |> tibble()
  table $ index <- rnames
  rm(result, rawtable, rnames)

  print(
    table[, c("index", "title", "species", "taxonomyid", "genome",
              "rdatadateadded", "sourcetype")]
  )
  stop()
}

taxoinfo <- list()
taxoinfo[["taxo"]] <- pargs $ taxo
taxoinfo[["s"]] <- pargs $ sgenes
taxoinfo[["g2m"]] <- pargs $ g2mgenes
saveRDS(taxoinfo, "taxo.rds")

if (file.exists(paste(gp_gencode, pargs $ taxo, "genomic.gff.table", sep = "/")) && # nolint
      file.exists(paste(gp_gencode, pargs $ taxo, "genomic.gff.exonlens", sep = "/"))) { # nolint

  if (file.exists("genome.rds")) file.remove("genome.rds")
  if (file.exists("exonlens.tsv")) file.remove("exonlens.tsv")
  if (!file.exists(paste(gp_gencode, pargs $ taxo, "genome.rds", sep = "/"))) {
    genes <- read.delim(
      paste(gp_gencode, pargs $ taxo, "genomic.gff.table", sep = "/"),
      sep = "\t", header = TRUE, comment.char = "#"
    )

    gene_table <- genes |> tibble()
    # the ENSEMBL gene ids have no duplicates.

    cat(blue(nrow(gene_table)), blue("genes in the annotation."), crlf)

    gene_table <- gene_table[!(genes $ gene_name == ""), ] # 48102 genes
    cat(blue(nrow(gene_table)), blue("genes with non-empty names."), crlf)

    # by now, non-duplicated columns are gene_id and gene_name.
    # we remove genes with duplicate id and assign chromosome identity

    chrfilter <- stringr::str_starts(gene_table $ .seqid, "chr")

    #. mitofilter <- stringr::str_starts(gene_table $ .seqid, pargs $ mt)
    #. chrfilter <- chrfilter & (!mitofilter)

    gene_table <- gene_table[chrfilter, ]
    cat(blue(nrow(gene_table)),
        blue("genes in primary chromosome assembly"), crlf)

    gene_table <- gene_table[!duplicated(gene_table $ gene_id), ]
    cat(blue(nrow(gene_table)),
        blue("genes with non-duplicated ensembl id"), crlf)

    # and here is the mitochondrial genome for mice in ensembl:

    # gene_id              .seqid .start  .end .strand gene_name gene_type
    # <chr>                <chr>   <int> <int> <chr>   <chr>     <chr>
    #  1 ENSMUSG00000064336.1 chrM        1    68 +       mt-Tf     Mt_tRNA
    #  2 ENSMUSG00000064337.1 chrM       70  1024 +       mt-Rnr1   Mt_rRNA
    #  3 ENSMUSG00000064338.1 chrM     1025  1093 +       mt-Tv     Mt_tRNA
    #  4 ENSMUSG00000064339.1 chrM     1094  2675 +       mt-Rnr2   Mt_rRNA
    #  5 ENSMUSG00000064340.1 chrM     2676  2750 +       mt-Tl1    Mt_tRNA
    #  6 ENSMUSG00000064341.1 chrM     2751  3707 +       mt-Nd1    protein_c...
    #  7 ENSMUSG00000064342.1 chrM     3706  3774 +       mt-Ti     Mt_tRNA
    #  8 ENSMUSG00000064343.1 chrM     3772  3842 -       mt-Tq     Mt_tRNA
    #  9 ENSMUSG00000064344.1 chrM     3845  3913 +       mt-Tm     Mt_tRNA
    # 10 ENSMUSG00000064345.1 chrM     3914  4951 +       mt-Nd2    protein_c...
    # 11 ENSMUSG00000064346.1 chrM     4950  5016 +       mt-Tw     Mt_tRNA
    # 12 ENSMUSG00000064347.1 chrM     5018  5086 -       mt-Ta     Mt_tRNA
    # 13 ENSMUSG00000064348.1 chrM     5089  5159 -       mt-Tn     Mt_tRNA
    # 14 ENSMUSG00000064349.1 chrM     5192  5257 -       mt-Tc     Mt_tRNA
    # 15 ENSMUSG00000064350.1 chrM     5260  5326 -       mt-Ty     Mt_tRNA
    # 16 ENSMUSG00000064351.1 chrM     5328  6872 +       mt-Co1    protein_c...
    # 17 ENSMUSG00000064352.1 chrM     6870  6938 -       mt-Ts1    Mt_tRNA
    # 18 ENSMUSG00000064353.1 chrM     6942  7011 +       mt-Td     Mt_tRNA
    # 19 ENSMUSG00000064354.1 chrM     7013  7696 +       mt-Co2    protein_c...
    # 20 ENSMUSG00000064355.1 chrM     7700  7764 +       mt-Tk     Mt_tRNA
    # 21 ENSMUSG00000064356.1 chrM     7766  7969 +       mt-Atp8   protein_c...
    # 22 ENSMUSG00000064357.1 chrM     7927  8607 +       mt-Atp6   protein_c...
    # 23 ENSMUSG00000064358.1 chrM     8607  9390 +       mt-Co3    protein_c...
    # 24 ENSMUSG00000064359.1 chrM     9391  9458 +       mt-Tg     Mt_tRNA
    # 25 ENSMUSG00000064360.1 chrM     9459  9806 +       mt-Nd3    protein_c...
    # 26 ENSMUSG00000064361.1 chrM     9808  9875 +       mt-Tr     Mt_tRNA
    # 27 ENSMUSG00000065947.1 chrM     9877 10173 +       mt-Nd4l   protein_c...
    # 28 ENSMUSG00000064363.1 chrM    10167 11544 +       mt-Nd4    protein_c...
    # 29 ENSMUSG00000064364.1 chrM    11546 11612 +       mt-Th     Mt_tRNA
    # 30 ENSMUSG00000064365.1 chrM    11613 11671 +       mt-Ts2    Mt_tRNA
    # 31 ENSMUSG00000064366.1 chrM    11671 11741 +       mt-Tl2    Mt_tRNA
    # 32 ENSMUSG00000064367.1 chrM    11742 13565 +       mt-Nd5    protein_c...
    # 33 ENSMUSG00000064368.1 chrM    13552 14070 -       mt-Nd6    protein_c...
    # 34 ENSMUSG00000064369.1 chrM    14071 14139 -       mt-Te     Mt_tRNA
    # 35 ENSMUSG00000064370.1 chrM    14145 15288 +       mt-Cytb   protein_c...
    # 36 ENSMUSG00000064371.1 chrM    15289 15355 +       mt-Tt     Mt_tRNA
    # 37 ENSMUSG00000064372.1 chrM    15356 15422 -       mt-Tp     Mt_tRNA

    mtids <- stringr::str_starts(gene_table $ .seqid, pargs $ mt)
    gene_table $ mito <- mtids
    cat(blue(sum(mtids)),
        blue("genes from mitochondrial genome"), crlf)

    # separate the ensembl accession and version
    gene_table <- gene_table |>
      tidyr::separate_wider_delim(
        gene_id, delim = ".",
        names = c("ensembl", "ensembl_version")
      )

    # convert to seurat-compatible gene names (_ is not allowed)
    gene_table $ seurat_names <-
      stringr::str_replace_all(gene_table $ gene_name, "_", "-")

    suppressPackageStartupMessages(
      {
        require(AnnotationHub)
        hub <- AnnotationHub::AnnotationHub()
        sql <- hub[[pargs $ annot]]
      }
    )

    go_columns <- columns(sql)
    go_columns <- go_columns[!(go_columns == "ENSEMBL")]

    all_genes <- rep("placeholder", length(go_columns))
    names(all_genes) <- c(go_columns)
    all_genes <- data.frame(all_genes) |> t()

    batch <- 5000
    ensembl <- gene_table $ ensembl

    for (x in seq(1, length(ensembl), batch)) {
      start <- x
      end <- x + batch - 1
      if (end > length(ensembl)) {
        end <- length(ensembl)
      }

      batch_list <- ensembl[start:end]
      batch_genes <- list()
      cat("extracting gene info from database orgdb: ", c(start, end), crlf)

      for (col in go_columns) {
        suppressWarnings(suppressMessages(
          query_list <- AnnotationDbi::select(
            sql, batch_list, col, verbose = FALSE, keytype = "ENSEMBL"
          )
        ))
        cols_data <- list()
        for (i in seq_along(batch_list)) {
          certain_col <- query_list[query_list $ ENSEMBL == batch_list[i], col]
          items_list <- certain_col |> unique() |> as.list()
          cols_data[[batch_list[i]]] <- do.call(paste, c(items_list, sep = ";"))
        }

        batch_genes <- cbind(batch_genes, cols_data)
      }

      all_genes <- rbind(all_genes, batch_genes)
    }

    row_name <- rownames(all_genes)
    all_genes <- all_genes |> as.data.table() |> tibble()
    all_genes <- all_genes[-1, ]
    all_genes $ ensembl <- row_name[-1]

    gene_table <- gene_table |> as.data.table() |> tibble()

    # by now the all_genes table has columns in lists. which will cause error
    # when ordering and deduplicating the rows, so we should unlist them to
    # string vectors for later operations

    for (cx in colnames(all_genes)) {
      all_genes[, cx] <- all_genes[, cx] |> unlist()
    }

    merged_table <- merge(
      gene_table, all_genes,
      all.x = TRUE, all.y = TRUE,
      by.x = "ensembl", by.y = "ensembl"
    )

    # this removal of duplicate need further investigation on why this occurs!
    merged_table <- merged_table[!duplicated(merged_table), ] |> tibble()
    cat(crlf)
    print(merged_table)

    newnames <- c()
    if (pargs $ rename) {
      cat(crlf)
      cat(stringr::str_wrap(paste(
        "make sure you have ", blue("entrez"), ", ", blue("refseq"), ", and ",
        blue("ensembl"), " in your data columns, since these are required for ",
        "gene assignment in the later steps.", sep = ""
      )), crlf)
      cat(blue("rename the gene table columns (<.> to keep original):"), crlf)

      for (cname in colnames(merged_table)) {
        cat("col", yellow(cname))
        newname <- read()
        if (newname == ".") newname <- cname
        newnames <- c(newnames, newname)
      }

      colnames(merged_table) <- newnames
    }

    saveRDS(
      merged_table,
      paste(gp_gencode, pargs $ taxo, "genome.rds", sep = "/")
    )
  }

  # create a link of the reference genome annotation R dataset to the current
  # working directory, as a sign that reference genome is specified.

  system(paste(
    "ln", "-s", paste(gp_gencode, pargs $ taxo, "genome.rds", sep = "/"),
    "genome.rds"
  ))

  system(paste(
    "ln", "-s",
    paste(gp_gencode, pargs $ taxo, "genomic.gff.exonlens", sep = "/"),
    "exonlens.tsv"
  ))

} else {
  cat(str_wrap(
    width = 80, indent = 0, exdent = 7, paste(
      red("error:"),
      "we can not find genomic.gff.table or genomic.gff.exonlens",
      "in the specified taxo directory, make sure you run gffindex,",
      "gffexonlen, and gffselect properly"
    )
  ))
}
