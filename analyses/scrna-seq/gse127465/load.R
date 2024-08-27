
library(Matrix)
library(SingleCellExperiment)

# read the read the file locations.
human_files = c(
  'geo/GSE127465/human/genes.tsv.gz',
  'geo/GSE127465/human/cells.tsv.gz',
  'geo/GSE127465/human/normalized.mtx.gz'
)

human_matrix <- Matrix::readMM(human_files[3])
human_gene_meta <- read.delim(human_files[1], comment.char = "#", header = FALSE)
human_cell_meta <- read.delim(human_files[2], comment.char = "#")

# this data has already been normalized and changed to log-2 based.

# the human dataset contains 54773 rows (cells) and 41861 columns (genes)
# however, we need to transpose that. rows should be genes and cols should be cells
human <- SingleCellExperiment(assays = list(
  normlog = human_matrix |> t()
))

dim(human)
rownames(human) <- human_gene_meta $ V1

colData(human) $ patient <- human_cell_meta $ Patient
colData(human) $ tissue <- human_cell_meta $ Tissue
colData(human) $ batch <- human_cell_meta $ Library
# this is the patient and tissue category combined
colData(human) $ sample <- human_cell_meta $ Barcoding.emulsion
colData(human) $ reads <- human_cell_meta $ Total.counts
colData(human) $ percent_mt <- human_cell_meta $ Percent.counts.from.mitochondrial.genes
colData(human) $ type <- human_cell_meta $ Most.likely.LM22.cell.type
colData(human) $ annot_maj <- human_cell_meta $ Major.cell.type
colData(human) $ annot_min <- human_cell_meta $ Minor.subset

# several pre-computed tags.
colData(human) $ is_blood <- human_cell_meta $ used_in_blood
colData(human) $ is_tumor <- human_cell_meta $ used_in_NSCLC_all_cells
colData(human) $ is_immune <- human_cell_meta $ used_in_NSCLC_and_blood_immune
colData(human) $ is_tumor_immune <- human_cell_meta $ used_in_NSCLC_immune
colData(human) $ is_tumor_not_immune <- human_cell_meta $ used_in_NSCLC_non_immune

human |> colData() |> colnames()

#  [1] "patient"             "tissue"              "batch"        "sample"
#  [5] "reads"               "percent_mt"          "type"         "annot_maj"
#  [9] "annot_min"           "is_blood"            "is_tumor"     "is_immune"
# [13] "is_tumor_immune"     "is_tumor_not_immune"

human $ type |> table()
human $ type |> length()
human $ annot_maj |> table()

colnames(human) <- paste(colData(human) $ sample, 1 : ncol(human), sep = ".")
saveRDS(human, "geo/GSE127465/sce.human.rds")

################################################################################

# the mouse dataset contains 28205 genes and 15939 cells.
mouse_files = c(
  'geo/GSE127465/mouse/genes.tsv.gz',
  'geo/GSE127465/mouse/cells.tsv.gz',
  'geo/GSE127465/mouse/normalized.mtx.gz'
)

mouse_matrix <- Matrix::readMM(mouse_files[3])
mouse_gene_meta <- read.delim(mouse_files[1], comment.char = "#", header = FALSE)
mouse_cell_meta <- read.delim(mouse_files[2], comment.char = "#")

mouse <- SingleCellExperiment(assays = list(
  normlog = mouse_matrix |> t()
))

rownames(mouse) <- mouse_gene_meta $ V1

colData(mouse) $ rep <- mouse_cell_meta $ Biological.replicate
colData(mouse) $ tissue <- mouse_cell_meta $ Tumor.or.healthy
colData(mouse) $ batch <- mouse_cell_meta $ Library
colData(mouse) $ reads <- mouse_cell_meta $ Total.counts
colData(mouse) $ percent_mt <- mouse_cell_meta $ Percent.counts.from.mitochondrial.genes
colData(mouse) $ type <- mouse_cell_meta $ Most.likely.Immgen.cell.type
colData(mouse) $ annot_maj <- mouse_cell_meta $ Major.cell.type
colData(mouse) $ annot_min <- mouse_cell_meta $ Minor.subset
colnames(mouse) <- paste(colData(mouse) $ batch, 1 : ncol(mouse), sep = ".")

saveRDS(mouse, "geo/GSE127465/sce.mouse.rds")
