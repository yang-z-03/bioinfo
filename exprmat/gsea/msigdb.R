
library(clusterProfiler)
library(data.table)
library(dplyr)
library(purrr)

# the pathways read here uses ncbi entrez id.
pathways <- read.gmt("msigdb/2023.2/pathways.gmt")
genes <- readRDS("refseq/9606/genes-ext.rds")
pathways <- merge(pathways, genes[, c("entrez", "refseq_gene")],
                  by.x = "gene", by.y = "entrez",
                  all.x = TRUE, all.y = FALSE) |>
  select(term, refseq_gene)

# all msigdb for human.
library(msigdbr)
hsa <- msigdbr(species = "Homo sapiens")
saveRDS(hsa |> as.data.frame(), "msigdb/2023.2/all.rds")

# another method to obtain the pathways subset of msigdb, rather than directly
# read the gmt files and convert names by yourself.
hsa <- readRDS("msigdb/2023.2/all.rds")
table(hsa $ gs_subcat)
pathways <- hsa |>
  filter(gs_subcat %in% c("CP:REACTOME")) |>
  select(gs_name, gene_symbol)

# prepare de genes (method 1) -------------------------------------------------
# using a concatenated table with classification
degenes <- readRDS("geo/GSE202374/merged-de.rds")
de <- degenes $ avg_log2FC
names(de) <- degenes $ gene
classification <- degenes $ cluster
rm(degenes)

# construct classification clustering
clusters <- map(0:30, function(x) {
  de[classification == x] |> sort(decreasing = TRUE)
})

cluster_names <- map(0:30, function(x) {
  de[classification == x] |> names()
})

# prepare de genes (method 2) -------------------------------------------------
# using diff-expr.R results
degenes <- readRDS("geo/GSE202374/merged-de.rds")
clusters <- map(0:30, function(x) {
  o <- degenes[, x + 1]
  names(o) <- rownames(degenes)
  o[abs(o) >= 1] |> sort(decreasing = TRUE)
})

cluster_names <- map(0:30, function(x) {
  o <- degenes[, x + 1]
  names(o) <- rownames(degenes)
  o[abs(o) >= 1] |> sort(decreasing = TRUE) |> names()
})

# running single sample enricher and gsea.
ID <- 3
opa <- enricher(clusters[[ID]] |> names(), TERM2GENE = pathways)
gse <- GSEA(clusters[[ID]], TERM2GENE = pathways, pvalueCutoff = 0.05)

library(enrichplot)
library(ggplot2)
library(ggraph)

dotplot(opa, showCategory = 30)

dotplot(gse, showCategory = 30)
upsetplot(gse)
cnetplot(gse, foldChange = clusters[[ID]], node_label = "gene")
heatplot(gse, foldChange = clusters[[ID]], showCategory = 5) + unify_theme_font()
ridgeplot(gse)
gseaplot(gse, geneSetID = 1, by = "runningScore")
gseaplot2(gse, geneSetID = 1)

gse2 <- pairwise_termsim(gse)
treeplot(gse2)
x <- emapplot(gse2, cex.params = list(category_label = 0)) +
  ggraph::geom_node_text(aes(label = c(
    'MHC Class II Antigen Presentation',
    'Adaptive Immune System',
    'Cytokine Signaling in Immune System',
    'Interferon Signaling'
  )), repel = TRUE, check_overlap = TRUE)
