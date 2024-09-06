
parser <- argparse::ArgumentParser(
  prog = "refer",
  description = "select an organism as reference genome"
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

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

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

if (file.exists(paste(gp_refseq, pargs $ taxo, "genomic.gff.table", sep = "/")) && # nolint
      file.exists(paste(gp_refseq, pargs $ taxo, "genomic.gff.exonlens", sep = "/"))) { # nolint

  if (file.exists("genome.rds")) file.remove("genome.rds")
  if (!file.exists(paste(gp_refseq, pargs $ taxo, "genome.rds", sep = "/"))) {
    genes <- read.delim(
      paste(gp_refseq, pargs $ taxo, "genomic.gff.table", sep = "/"),
      sep = "\t", header = TRUE, comment.char = "#"
    )

    id_pattern <- "^(.*;)?( *GeneID: *([^,]+))(.*)$"
    selected_id_lines <- grep(id_pattern, genes $ Dbxref,
                              perl = TRUE, ignore.case = TRUE)

    ncbi_id <- sub(id_pattern, "\\3",
                   genes[selected_id_lines, "Dbxref"],
                   perl = TRUE, ignore.case = TRUE)

    genes <- genes[selected_id_lines, ]
    gene_table <- cbind(genes, ncbi_id)

    suppressPackageStartupMessages(
      {
        require(AnnotationHub)
        hub <- AnnotationHub::AnnotationHub()
        sql <- hub[[pargs $ annot]]
      }
    )
    go_columns <- columns(sql)
    go_columns <- go_columns[!(go_columns == "ENTREZID")]

    all_genes <- rep("placeholder", length(go_columns))
    names(all_genes) <- c(go_columns)
    all_genes <- data.frame(all_genes) |> t()

    batch <- 5000
    entrez <- ncbi_id

    for (x in seq(1, length(entrez), batch)) {
      start <- x
      end <- x + batch - 1
      if (end > length(entrez)) {
        end <- length(entrez)
      }

      batch_list <- entrez[start:end]
      batch_genes <- list()
      cat("extracting gene info from database orgdb: ", c(start, end), crlf)

      for (col in go_columns) {
        suppressMessages(
          query_list <- select(sql, batch_list, col, verbose = FALSE)
        )
        cols_data <- list()
        for (i in seq_along(batch_list)) {
          certain_col <- query_list[query_list $ ENTREZID == batch_list[i], col]
          items_list <- certain_col |> unique() |> as.list()
          cols_data[[batch_list[i]]] <- do.call(paste, c(items_list, sep = ";"))
        }

        batch_genes <- cbind(batch_genes, cols_data)
      }

      all_genes <- rbind(all_genes, batch_genes)
    }

    all_genes <- all_genes[-1, ]
    row_name <- rownames(all_genes)
    all_genes <- all_genes |> as.data.table() |> tibble()
    all_genes $ ncbi_id <- row_name

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
      by.x = "ncbi_id", by.y = "ncbi_id"
    )

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
      paste(gp_refseq, pargs $ taxo, "genome.rds", sep = "/")
    )
  }

  # create a link of the reference genome annotation R dataset to the current
  # working directory, as a sign that reference genome is specified.

  system(paste(
    "ln", "-s", paste(gp_refseq, pargs $ taxo, "genome.rds", sep = "/"),
    "genome.rds"
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
