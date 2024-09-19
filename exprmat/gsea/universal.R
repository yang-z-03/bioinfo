
library(clusterProfiler)
library(data.table)
library(dplyr)

degenes <- readRDS("geo/GSE202374/merged-differential-expressions.rds")
genes <- readRDS("refseq/9606/genes-ext.rds")

# here, in degenes, we have the log2-fold-change value, and the gene name.
# we can obtain the go-all and kegg fields from genes to manually construct
# a gene set definition.

# prepare the fold changes data.

de <- degenes $ avg_log2FC
names(de) <- degenes $ gene
de <- sort(de, decreasing = TRUE)
rm(degenes)

# prepare kegg gene sets.

nrow(genes)

construct_kegg <- function(tbgene) {
  expand_names <- list()
  expand_kegg <- list()
  for(line in 1 : nrow(tbgene)) {
    gene_name <- tbgene[line, ] $ refseq_name
    gene_kegg <- tbgene[line, ] $ kegg

    if (gene_kegg == "NA") {}
    else {
      expand_kegg[[line]] <- strsplit(gene_kegg, ";") |> unlist()
      expand_names[[line]] <- rep(gene_name, length(expand_kegg[[line]]))
    }
  }

  data.frame(list(
    gene = expand_names |> unlist(),
    kegg = expand_kegg |> unlist()
  ))
}

construct_go <- function(tbgene) {
  expand_names <- list()
  expand_go <- list()
  for(line in 1 : nrow(tbgene)) {
    gene_name <- tbgene[line, ] $ refseq_name
    gene_go <- tbgene[line, ] $ go_all

    if (gene_go == "NA") {}
    else {
      expand_go[[line]] <- strsplit(gene_go, ";") |> unlist()
      expand_names[[line]] <- rep(gene_name, length(expand_go[[line]]))
    }
  }

  data.frame(list(
    gene = expand_names |> unlist(),
    go = expand_go |> unlist()
  ))
}

kegg_set <- construct_kegg(genes)
go_set <- construct_go(genes)

# over-presentation kegg
opa <- enricher(de |> names(), TERM2GENE = kegg_set |> dplyr::select(kegg, gene))
gse <- GSEA(de, TERM2GENE = kegg_set |> dplyr::select(kegg, gene))

# similarly, you can do opa and gse for go.
opa <- enricher(de |> names(), TERM2GENE = go_set |> dplyr::select(go, gene))
gse <- GSEA(de, TERM2GENE = go_set |> dplyr::select(go, gene))

library(enrichplot)

dotplot(gse, showCategory = 30)
upsetplot(gse)
ridgeplot(gse)
gseaplot(gse, geneSetID = 1, by = "runningScore")
gseaplot2(gse, geneSetID = 1:2)
