
library(clusterProfiler)

genes <- readRDS("refseq/9606/genes.rds") |> data.frame()

id_pattern <- "^(.*;)?( *GeneID: *([^,]+))(.*)$"
selected_id_lines <- grep(id_pattern, genes $ xref,
                          perl = TRUE, ignore.case = TRUE)

ncbi_id <- sub(id_pattern, "\\3",
               genes[selected_id_lines, "xref"],
               perl = TRUE, ignore.case = TRUE)

hgnc_pattern <- "^(.*,)?( *HGNC: *([^,]+))(.*)$"
selected_hgnc_lines <- grep(hgnc_pattern, genes $ xref,
                            perl = TRUE, ignore.case = TRUE)

hgnc_id <- sub(hgnc_pattern, "\\3",
               genes[selected_hgnc_lines, "xref"],
               perl = TRUE, ignore.case = TRUE)

kegg <- clusterProfiler::bitr_kegg(
  ncbi_id, fromType = "ncbi-geneid",
  toType = "kegg",
  organism = "hsa"
)

ncbi_protid <- clusterProfiler::bitr_kegg(
  ncbi_id, fromType = "ncbi-geneid",
  toType = "ncbi-proteinid",
  organism = "hsa"
)

ncbi_uniprot <- clusterProfiler::bitr_kegg(
  ncbi_id, fromType = "ncbi-geneid",
  toType = "uniprot",
  organism = "hsa"
)

merged_kegg <- kegg |>
  merge(ncbi_protid, all = TRUE) |>
  merge(ncbi_uniprot, all = TRUE)

library(AnnotationHub)
hub <- AnnotationHub()

query(hub, c("sapiens", "orgdb"))

# AnnotationHub with 1 record
# # snapshotDate(): 2024-04-30
# # names(): AH116710
# # $dataprovider: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/
# # $species: Homo sapiens
# # $rdataclass: OrgDb
# # $rdatadateadded: 2024-04-02
# # $title: org.Hs.eg.db.sqlite
# # $description: NCBI gene ID based annotations about Homo sapiens
# # $taxonomyid: 9606
# # $genome: NCBI genomes
# # $sourcetype: NCBI/ensembl
# # $sourceurl: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/, ftp://ftp.ensembl.org/...
# # $sourcesize: NA
# # $tags: c("NCBI", "Gene", "Annotation")
# # retrieve record with 'object[["AH116710"]]'

hsa <- hub[["AH116710"]]
columns(hsa)

# ncbi gene id is the same thing as entrez id.
# construct the gene list with entrez id included.

gene_table <- cbind(genes, ncbi_id)
head(gene_table)
colnames(gene_table) <- c(
  "refseq-seqname", "refseq-source", "refseq-feature", "refseq-start",
  "refseq-end", "refseq-score", "refseq-strand", "refseq-frame",
  "refseq-id", "refseq-xref", "refseq-name", "refseq-description",
  "refseq-genbank", "refseq-gene", "refseq-biotype", "entrez-id"
)

# entrez id is the major key.

go_columns <- c(
  "ACCNUM", "ALIAS", "ENSEMBL", "ENSEMBLPROT", "ENSEMBLTRANS",
  "ENZYME", "EVIDENCE", "EVIDENCEALL", "GENENAME",  "GENETYPE", "GO", "GOALL",
  "IPI", "MAP", "OMIM", "ONTOLOGY",  "ONTOLOGYALL", "PATH", "PFAM", "PMID",
  "PROSITE", "REFSEQ", "SYMBOL",  "UCSCKG", "UNIPROT"
)

select(
  hsa, c("1"), c("SYMBOL")
)

all_genes <- rep("placeholder", length(go_columns))
names(all_genes) <- c(go_columns)
all_genes <- data.frame(all_genes) |> t()

batch <- 10000
entrez <- gene_table $ `entrez-id`

for (x in seq(1, length(entrez), batch)) {
  start <- x
  end <- x + batch - 1
  if (end > length(entrez)) {
    end <- length(entrez)
  }

  batch_list <- entrez[start:end]
  batch_genes <- list()
  print(c(start, end))

  for (col in go_columns) {
    query_list <- select(hsa, batch_list, col, verbose = FALSE)
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
library(data.table)
all_genes <- all_genes |> as.data.table()
all_genes $ ncbi_id <- row_name

# merge all_genes columns with gene_table.
# gene_table (48104), all_genes(45593)

gene_table <- gene_table |> as.data.table()

merged_table <- merge(
  gene_table, all_genes,
  all.x = TRUE, all.y = TRUE,
  by.x = "entrez-id", by.y = "ncbi_id"
) # 56254

gene_table $ `entrez-id` |> duplicated() |> table()

# FALSE  TRUE
# 42651  5453

merged_table $ `entrez-id` |> duplicated() |> table()

# FALSE  TRUE
# 42651 13603

merged_table[merged_table $ `entrez-id` |> duplicated()] $ `entrez-id` |> head()

# [1] "10000"     "10000"     "10000"     "100009603" "100009603" "100009603"

merged_table[`entrez-id` == "100009603"]

gene_table |> duplicated() |> table()

colnames(merged_table) <- c(
  "entrez",

  "refseq_seqname", "refseq_source", "refseq_feature", "refseq_start",
  "refseq_end", "refseq_score", "refseq_strand", "refseq_frame",
  "refseq_id", "refseq_xref", "refseq_name", "refseq_description",
  "refseq_genbank", "refseq_gene", "refseq_biotype",

  "genbank_accession", "symbol_alias",
  "ensembl", "ensembl_prot", "ensembl_transcript",
  "ec",
  "go_evidence", "go_evidence_all", "name", "type", "go", "go_all",
  "ipi_accession", "map_location",
  "omim", "go_ontology",  "go_ontology_all",
  "kegg_pathway", "pfam", "pmid",
  "prosite", "refseq", "symbol",  "ucsc_known_gene", "uniprot"
)

cast <- merged_table[, genbank_accession := as.character(genbank_accession)]
cast <- cast[, symbol_alias := as.character(symbol_alias)]
cast <- cast[, ensembl := as.character(ensembl)]
cast <- cast[, ensembl_prot := as.character(ensembl_prot)]
cast <- cast[, ensembl_transcript := as.character(ensembl_transcript)]
cast <- cast[, ec := as.factor(as.character(ec))]
cast <- cast[, go_evidence := as.character(go_evidence)]
cast <- cast[, go_evidence_all := as.character(go_evidence_all)]
cast <- cast[, name := as.character(name)]
cast <- cast[, type := as.factor(as.character(type))]
cast <- cast[, go := as.character(go)]
cast <- cast[, go_all := as.character(go_all)]
cast <- cast[, ipi_accession := as.character(ipi_accession)]
cast <- cast[, map_location := as.character(map_location)]
cast <- cast[, omim := as.character(omim)]
cast <- cast[, go_ontology := as.character(go_ontology)]
cast <- cast[, go_ontology_all := as.character(go_ontology_all)]
cast <- cast[, kegg_pathway := as.character(kegg_pathway)]
cast <- cast[, pfam := as.character(pfam)]
cast <- cast[, pmid := as.character(pmid)]
cast <- cast[, prosite := as.character(prosite)]
cast <- cast[, refseq := as.character(refseq)]
cast <- cast[, symbol := as.character(symbol)]
cast <- cast[, ucsc_known_gene := as.character(ucsc_known_gene)]
cast <- cast[, uniprot := as.character(uniprot)]

cast |> duplicated() |> table()
cast <- cast[!duplicated(cast), ]

saveRDS(cast, file = "refseq/9606/genes-ext.rds")

# > colnames(cast)
#
#  [1] "entrez"             "refseq_seqname"     "refseq_source"
#  [4] "refseq_feature"     "refseq_start"       "refseq_end"
#  [7] "refseq_score"       "refseq_strand"      "refseq_frame"
# [10] "refseq_id"          "refseq_xref"        "refseq_name"
# [13] "refseq_description" "refseq_genbank"     "refseq_gene"
# [16] "refseq_biotype"     "genbank_accession"  "symbol_alias"
# [19] "ensembl"            "ensembl_prot"       "ensembl_transcript"
# [22] "ec"                 "go_evidence"        "go_evidence_all"
# [25] "name"               "type"               "go"
# [28] "go_all"             "ipi_accession"      "map_location"
# [31] "omim"               "go_ontology"        "go_ontology_all"
# [34] "kegg_pathway"       "pfam"               "pmid"
# [37] "prosite"            "refseq"             "symbol"
# [40] "ucsc_known_gene"    "uniprot"