
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

parser $ add_argument(
  "-m", dest = "mt", type = "character",
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

if (file.exists(paste(gp_refseq, pargs $ taxo, "genomic.gff.table", sep = "/")) && # nolint
      file.exists(paste(gp_refseq, pargs $ taxo, "genomic.gff.exonlens", sep = "/"))) { # nolint

  if (file.exists("genome.rds")) file.remove("genome.rds")
  if (file.exists("exonlens.tsv")) file.remove("exonlens.tsv")
  if (!file.exists(paste(gp_refseq, pargs $ taxo, "genome.rds", sep = "/"))) {
    genes <- read.delim(
      paste(gp_refseq, pargs $ taxo, "genomic.gff.table", sep = "/"),
      sep = "\t", header = TRUE, comment.char = "#"
    )

    id_pattern <- "^(.*;)?( *GeneID: *([^,]+))(.*)$"
    selected_id_lines <- grep(id_pattern, genes $ Dbxref,
                              perl = TRUE, ignore.case = TRUE)

    # commonly, every gene should have this field. and this selection will
    # not be effective normally.

    ncbi_id <- sub(id_pattern, "\\3",
                   genes[selected_id_lines, "Dbxref"],
                   perl = TRUE, ignore.case = TRUE)

    # here, there may be some notice points for the ENTREZ annotations:
    #
    # (1) 4 genes have misteriously no symbols and id names. however with a
    #     .type == gene specification. maybe some mistakes.
    #
    # (2) refseq column 'gene' and 'name' are completely identical, both as
    #     the standard hugo gene names.

    genes <- genes[selected_id_lines, ]
    gene_table <- cbind(genes, ncbi_id) # 48106 genes
    cat(blue(nrow(gene_table)), blue("genes in the annotation."), crlf)

    # remove those with empty symbol names and IDs (just 4, throw away)
    gene_table <- gene_table[!(genes $ ID == ""), ] # 48102 genes
    cat(blue(nrow(gene_table)), blue("genes with non-empty names."), crlf)

    # by now, the only non-duplicated column is ID.

    # it seems that the ncbi refseq is the least duplicated system that can be
    # used though. in concerns of transcriptome mapping, we may just need the
    # gene annotations on the well-defined chromosomes. by confining the .seqid
    # field to those starting with NC_, we can obtain the 24 chromosomes (22XY)
    # of human for example, and a special chromosome of mitochondria.

    # for human, the mt genome is marked as: "NC_012920.1"

    # after the filtering step for major chromosomes, we still have 39 gene
    # name duplications listed below (for human): (and only 34 duplicated
    # ncbi entrez ids (since the TRNAA/V are only five))

    #    ncbi_id   .seqid       .type    .start      .end .strand gene
    #    <chr>     <chr>        <chr>     <int>     <int> <chr>   <chr>
    #  1 107985615 NC_000001.11 gene  145161108 145161178 +       TRNAV-CAC
    #  2 124901561 NC_000006.12 gene   26771080  26771152 -       TRNAA-AGC
    #  3 124901562 NC_000006.12 gene   26796209  26796281 -       TRNAA-AGC
    #  4 124901563 NC_000006.12 gene   26814339  26814411 -       TRNAA-AGC
    #  5 124901564 NC_000006.12 gene   26819109  26819181 -       TRNAA-AGC
    #  6 55344     NC_000024.10 gene     276356    303356 +       PLCXD1
    #  7 8225      NC_000024.10 gene     304759    318796 -       GTPBP6
    #  8 283981    NC_000024.10 gene     319145    321332 +       LINC00685
    #  9 28227     NC_000024.10 gene     333933    386907 -       PPP2R3B
    # 10 102724521 NC_000024.10 gene     430021    472756 +       LOC102724521
    # 11 6473      NC_000024.10 gene     624344    659411 +       SHOX
    # 12 64109     NC_000024.10 gene    1190490   1212649 -       CRLF2
    # 13 1438      NC_000024.10 gene    1268814   1325218 +       CSF2RA
    # 14 124905238 NC_000024.10 gene    1293546   1293998 -       LOC124905238
    # 15 100500894 NC_000024.10 gene    1293918   1293992 +       MIR3690
    # 16 3563      NC_000024.10 gene    1336785   1382689 +       IL3RA
    # 17 101928032 NC_000024.10 gene    1336972   1378476 -       LOC101928032
    # 18 293       NC_000024.10 gene    1386152   1392113 -       SLC25A6
    # 19 124900597 NC_000024.10 gene    1391993   1396881 +       LOC124900597
    # 20 751580    NC_000024.10 gene    1397025   1399412 +       LINC00106
    # 21 80161     NC_000024.10 gene    1400531   1415421 +       ASMTL-AS1
    # 22 8623      NC_000024.10 gene    1403139   1453756 -       ASMTL
    # 23 286530    NC_000024.10 gene    1462581   1537185 -       P2RY8
    # 24 8227      NC_000024.10 gene    1591604   1602520 +       AKAP17A
    # 25 438       NC_000024.10 gene    1615059   1643081 +       ASMT
    # 26 105373105 NC_000024.10 gene    1732419   1755984 +       LINC02968
    # 27 105379413 NC_000024.10 gene    1771053   1784602 +       LOC105379413
    # 28 107985677 NC_000024.10 gene    1879515   1885827 -       LOC107985677
    # 29 207063    NC_000024.10 gene    2219506   2500976 -       DHRSX
    # 30 124905239 NC_000024.10 gene    2320639   2337671 +       LOC124905239
    # 31 9189      NC_000024.10 gene    2486435   2500976 -       ZBED1
    # 32 101928092 NC_000024.10 gene    2566029   2609240 -       LINC03112
    # 33 102464837 NC_000024.10 gene    2609191   2609254 +       MIR6089
    # 34 100359394 NC_000024.10 gene    2612991   2615347 -       LINC00102
    # 35 4267      NC_000024.10 gene    2691295   2741309 +       CD99
    # 36 10251     NC_000024.10 gene   56923423  56968979 +       SPRY3
    # 37 6845      NC_000024.10 gene   57067865  57130289 +       VAMP7
    # 38 3581      NC_000024.10 gene   57184072  57197869 +       IL9R
    # 39 100128260 NC_000024.10 gene   57201084  57203350 -       WASIR1

    # the reason of duplication is fairly clear for human:
    #
    # (1) the genes encoding tRNAs just have multiple positions on the
    #     chromosome, these genes have different entrez ids in the ncbi database
    #     but just must have the same name. the behavior on dealing with such
    #     genes in transcriptomic sequencing will be merging the genes with
    #     the same name and add their expression together.
    #
    # (2) the genes are alleles between the X and Y chromosome. this is the
    #     only included paired chromosome in the genome. so they must have same
    #     genes. they have the same ncbi entrez numbers since they are alleles.

    chrfilter <- stringr::str_starts(gene_table $ .seqid, "NC")
    gene_table <- gene_table[chrfilter, ]
    cat(blue(nrow(gene_table)), blue("genes in primary chromosome assembly"), crlf)

    gene_table <- gene_table[!duplicated(gene_table $ ncbi_id), ]
    cat(blue(nrow(gene_table)), blue("genes with non-duplicated entrez id"), crlf)

    # and here is the mitochondrial genome for human:

    #    ncbi_id .seqid      .start  .end .strand gene  gene_biotype
    #    <chr>   <chr>        <int> <int> <chr>   <chr> <chr>
    #  1 4558    NC_012920.1    577   647 +       TRNF  tRNA
    #  2 4549    NC_012920.1    648  1601 +       RNR1  rRNA
    #  3 4577    NC_012920.1   1602  1670 +       TRNV  tRNA
    #  4 4550    NC_012920.1   1671  3229 +       RNR2  rRNA
    #  5 4567    NC_012920.1   3230  3304 +       TRNL1 tRNA
    #  6 4535    NC_012920.1   3307  4262 +       ND1   protein_coding
    #  7 4565    NC_012920.1   4263  4331 +       TRNI  tRNA
    #  8 4572    NC_012920.1   4329  4400 -       TRNQ  tRNA
    #  9 4569    NC_012920.1   4402  4469 +       TRNM  tRNA
    # 10 4536    NC_012920.1   4470  5511 +       ND2   protein_coding
    # 11 4578    NC_012920.1   5512  5579 +       TRNW  tRNA
    # 12 4553    NC_012920.1   5587  5655 -       TRNA  tRNA
    # 13 4570    NC_012920.1   5657  5729 -       TRNN  tRNA
    # 14 4511    NC_012920.1   5761  5826 -       TRNC  tRNA
    # 15 4579    NC_012920.1   5826  5891 -       TRNY  tRNA
    # 16 4512    NC_012920.1   5904  7445 +       COX1  protein_coding
    # 17 4574    NC_012920.1   7446  7514 -       TRNS1 tRNA
    # 18 4555    NC_012920.1   7518  7585 +       TRND  tRNA
    # 19 4513    NC_012920.1   7586  8269 +       COX2  protein_coding
    # 20 4566    NC_012920.1   8295  8364 +       TRNK  tRNA
    # 21 4509    NC_012920.1   8366  8572 +       ATP8  protein_coding
    # 22 4508    NC_012920.1   8527  9207 +       ATP6  protein_coding
    # 23 4514    NC_012920.1   9207  9990 +       COX3  protein_coding
    # 24 4563    NC_012920.1   9991 10058 +       TRNG  tRNA
    # 25 4537    NC_012920.1  10059 10404 +       ND3   protein_coding
    # 26 4573    NC_012920.1  10405 10469 +       TRNR  tRNA
    # 27 4539    NC_012920.1  10470 10766 +       ND4L  protein_coding
    # 28 4538    NC_012920.1  10760 12137 +       ND4   protein_coding
    # 29 4564    NC_012920.1  12138 12206 +       TRNH  tRNA
    # 30 4575    NC_012920.1  12207 12265 +       TRNS2 tRNA
    # 31 4568    NC_012920.1  12266 12336 +       TRNL2 tRNA
    # 32 4540    NC_012920.1  12337 14148 +       ND5   protein_coding
    # 33 4541    NC_012920.1  14149 14673 -       ND6   protein_coding
    # 34 4556    NC_012920.1  14674 14742 -       TRNE  tRNA
    # 35 4519    NC_012920.1  14747 15887 +       CYTB  protein_coding
    # 36 4576    NC_012920.1  15888 15953 +       TRNT  tRNA
    # 37 4571    NC_012920.1  15956 16023 -       TRNP  tRNA

    mtids <- stringr::str_starts(gene_table $ .seqid, pargs $ mt)
    is_mito <- rep(FALSE, length(gene_table $ .seqid))
    is_mito[mtids] <- TRUE
    gene_table $ mito <- is_mito

    # convert to seurat-compatible gene names (_ is not allowed)
    gene_table $ seurat_names <-
      stringr::str_replace_all(gene_table $ gene, "_", "-")

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
    entrez <- gene_table $ ncbi_id

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
        suppressWarnings(suppressMessages(
          query_list <- AnnotationDbi::select(
            sql, batch_list, col, verbose = FALSE
          )
        ))
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

    row_name <- rownames(all_genes)
    all_genes <- all_genes |> as.data.table() |> tibble()
    all_genes <- all_genes[-1, ]
    all_genes $ ncbi_id <- row_name[-1]

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
      paste(gp_refseq, pargs $ taxo, "genome.rds", sep = "/")
    )
  }

  # create a link of the reference genome annotation R dataset to the current
  # working directory, as a sign that reference genome is specified.

  system(paste(
    "ln", "-s", paste(gp_refseq, pargs $ taxo, "genome.rds", sep = "/"),
    "genome.rds"
  ))

  system(paste(
    "ln", "-s",
    paste(gp_refseq, pargs $ taxo, "genomic.gff.exonlens", sep = "/"),
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
