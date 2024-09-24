
parser <- argparse::ArgumentParser(
  prog = "qc",
  description = "generate quality control metrics and perform filtered result"
)

parser $ add_argument(
  "-m", dest = "mtpct", type = "double", default = -1,
  help = "threshold for high mitochondrial gene percentage"
)

parser $ add_argument(
  "-r", dest = "rbpct", type = "double", default = -1,
  help = "threshold for low ribosomal gene percentage"
)

parser $ add_argument(
  "-l", dest = "ldepth", type = "integer", default = -1,
  help = "threshold for low read counts per cell"
)

parser $ add_argument(
  "-u", dest = "udepth", type = "integer", default = -1,
  help = "threshold for excess read counts per cell (presumbly doublets)"
)

parser $ add_argument(
  "-g", dest = "expgene", type = "integer", default = 0,
  help = "lower threshold for minimum gene detection"
)

parser $ add_argument(
  "-d", dest = "scrublet", action = "store_true", default = FALSE,
  help = "run scrublet estimation of doublets"
)

parser $ add_argument(
  "-e", dest = "expr", type = "double", default = 0.1,
  help = "minimum expressing cell percentage per gene"
)

parser $ add_argument(
  "-t", dest = "minexpr", type = "integer", default = 3,
  help = "minimum expressing transcript counting as expression"
)

parser $ add_argument(
  "-p", dest = "plot", action = "store_true", default = FALSE,
  help = "save plottings"
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

if (!shared[["is_loaded"]]) {
  cat("you should run 'read' or 'read10x' first to generate the uniform",
      "expression matrix.", crlf)
  stop()
}

if (!dir.exists("qc"))
  dir.create("qc")

# the most common quality control is to filter out
#
# 1.  cells with too few genes detected. they usually represent cells which are
#     not sequenced deep enough for reliable characterization.
#
# 2.  cells with too many genes detected. they may represent doublets or
#     multiplets (i.e. two or more cells in the same droplet, therefore sharing
#     the same cell barcode).
#
# 3.  cells with high mitochondrial transcript percentage. as most of the
#     scrna-seq experiments use oligo-t to capture mrnas, mitochondrial
#     transcripts should be relatively under-representative due to their lack
#     of poly-a tails, but it is unavoidable that some mitochondrial transcripts
#     are captured. meanwhile, there is also some evidence that stable poly-a
#     tails exist in some mitochondrial transcripts but serve as a marker for
#     degradation. together, cells with high mitochondrial transcript percentage
#     likely represent cells under stress (e.g. hypoxia) which produce a lot of
#     mitochondria, or which produce an abnormally high amount of truncated
#     mitochondrial transcripts.

expr_count <- readRDS("features/matrix.rds")
sample_meta <- readRDS("features/samples-meta.rds")
genes_meta <- readRDS("features/genes-meta.rds")
taxo <- readRDS("taxo.rds")

sample_meta $ reads <- colSums(expr_count)

expr_mt <- expr_count[genes_meta $ mito, ]
expr_rb <- expr_count[grep("^RP[LS]", genes_meta $ gene, ignore.case = TRUE), ]

sample_meta $ pct_ribo <- colSums(expr_rb) / pull(sample_meta, "reads")
sample_meta $ pct_mito <- colSums(expr_mt) / pull(sample_meta, "reads")

sample_meta $ detection <- colSums(expr_count > 0)

suppressPackageStartupMessages(
  source(paste(gp_base, "scrublet.R", sep = "/"))
)

if (pargs $ scrublet) {
  d <- scrublet_matrix(NULL, make_matrix(expr_count),
                       return_results_only = TRUE)
  sample_meta $ doublet_scores <- d[["doublet_scores"]]
  sample_meta $ is_doublet <- d[["predicted_doublets"]]
}

upper_mito <- pargs $ mtpct
if (upper_mito < 0)
  upper_mito <- (isOutlier(sample_meta $ pct_mito, type = "higher") |>
                   attr("thresholds"))["higher"]

lower_ribo <- pargs $ rbpct

lower_depth <- pargs $ ldepth
if (lower_depth < 0)
  lower_depth <- (isOutlier(sample_meta $ reads, type = "both") |>
                    attr("thresholds"))["lower"]

upper_depth <- pargs $ udepth
if (upper_depth < 0)
  upper_depth <- (isOutlier(sample_meta $ reads, type = "both") |>
                    attr("thresholds"))["higher"]

sample_meta $ qc <- (
  sample_meta $ pct_mito <= upper_mito &
    sample_meta $ reads >= lower_depth &
    sample_meta $ reads <= upper_depth &
    sample_meta $ detection >= pargs $ expgene
)

if (lower_ribo > 0)
  sample_meta $ qc <- sample_meta $ qc & (sample_meta $ pct_ribo >= lower_ribo)

if (pargs $ scrublet) {
  sample_meta $ qc <- sample_meta $ qc & (!sample_meta $ is_doublet)
  cat(crlf)
}

cat(green("per-sample qc:"), crlf)
print(table(sample_meta $ qc))
cat(crlf)


# filtering low-expression genes.

expr_cells <- rowSums(expr_count > pargs $ minexpr)
expr_thresh <- as.integer(ncol(expr_count) * pargs $ expr / 100)
rowfilter <- expr_cells >= expr_thresh

cat(green("per-gene qc:"), crlf)
print(table(expr_cells >= expr_thresh))
cat(crlf)

if (pargs $ plot) suppressMessages({

  if (lower_ribo < 0) lower_ribo <- 0

  h1 <- ggplot2::ggplot(sample_meta, aes(x = pct_ribo)) +
    ggplot2::geom_histogram() +
    geom_vline(aes(xintercept = lower_ribo),
               color = "#bd0000", linetype = "dashed", size = 0.5) +
    ggplot2::labs(x = "Ribosomal percentage", y = "Cell counts") +
    ggplot2::scale_x_continuous(expand = c(0.05, 0.05)) +
    unify_theme_font(base_family = "Myriad Pro Cond")

  h2 <- ggplot2::ggplot(sample_meta, aes(x = pct_mito)) +
    ggplot2::geom_histogram() +
    geom_vline(aes(xintercept = upper_mito),
               color = "#bd0000", linetype = "dashed", size = 0.5) +
    ggplot2::labs(x = "Mitochondrial percentage", y = "Cell counts") +
    ggplot2::scale_x_continuous(expand = c(0.05, 0.05)) +
    unify_theme_font(base_family = "Myriad Pro Cond")

  h3 <- ggplot2::ggplot(sample_meta, aes(x = reads)) +
    ggplot2::geom_histogram() +
    geom_vline(aes(xintercept = upper_depth),
               color = "#bd0000", linetype = "dashed", size = 0.5) +
    geom_vline(aes(xintercept = lower_depth),
               color = "#bd0000", linetype = "dashed", size = 0.5) +
    ggplot2::labs(x = "Total reads (sequencing depth)", y = "Cell counts") +
    ggplot2::scale_x_continuous(expand = c(0.05, 0.05)) +
    unify_theme_font(base_family = "Myriad Pro Cond")

  h4 <- ggplot2::ggplot(data.frame(.c = expr_cells), aes(x = .c)) +
    ggplot2::geom_histogram() +
    geom_vline(aes(xintercept = expr_thresh),
               color = "#bd0000", linetype = "dashed", size = 0.5) +
    ggplot2::labs(x = "Expressing cells", y = "Genes") +
    ggplot2::scale_x_continuous(expand = c(0.05, 0.05)) +
    unify_theme_font(base_family = "Myriad Pro Cond")

  s1 <- ggplot2::ggplot(sample_meta, aes(x = reads, y = pct_mito)) +
    ggplot2::geom_point(alpha = 0.5) +
    geom_vline(aes(xintercept = lower_depth),
               color = "#bd0000", linetype = "dashed", linewidth = 0.5) +
    geom_vline(aes(xintercept = upper_depth),
               color = "#bd0000", linetype = "dashed", linewidth = 0.5) +
    geom_hline(aes(yintercept = upper_mito),
               color = "#bd0000", linetype = "dashed", linewidth = 0.5) +
    ggplot2::labs(x = "Total reads (sequencing depth)",
                  y = "Mitochondrial percentage") +
    ggplot2::scale_x_continuous(expand = c(0.05, 0.05)) +
    unify_theme_font(base_family = "Myriad Pro Cond")

  s2 <- ggplot2::ggplot(sample_meta, aes(x = reads, y = detection)) +
    ggplot2::geom_point(alpha = 0.5) +
    geom_vline(aes(xintercept = lower_depth),
               color = "#bd0000", linetype = "dashed", linewidth = 0.5) +
    geom_vline(aes(xintercept = upper_depth),
               color = "#bd0000", linetype = "dashed", linewidth = 0.5) +
    ggplot2::labs(x = "Total reads (sequencing depth)",
                  y = "Gene detections") +
    ggplot2::scale_x_continuous(expand = c(0.05, 0.05)) +
    unify_theme_font(base_family = "Myriad Pro Cond")

  ggsave("qc/pct-ribo.png", h1, units = "in",
         width = 5, height = 2.5, dpi = 600)
  ggsave("qc/pct-mito.png", h2, units = "in",
         width = 5, height = 2.5, dpi = 600)
  ggsave("qc/reads.png", h3, units = "in",
         width = 5, height = 2.5, dpi = 600)
  ggsave("qc/genes.png", h4, units = "in",
         width = 5, height = 2.5, dpi = 600)

  ggsave("qc/reads-vs-mito.png", s1, units = "in",
         width = 5, height = 5, dpi = 600)
  ggsave("qc/reads-vs-detections.png", s2, units = "in",
         width = 5, height = 5, dpi = 600)

})

expr_count <- expr_count[, sample_meta $ qc]
sample_meta <- sample_meta[sample_meta $ qc, ]

expr_count <- expr_count[rowfilter, ]
genes_meta <- genes_meta[rowfilter, ]

saveRDS(sample_meta, "qc/samples-meta.rds")
saveRDS(genes_meta, "qc/genes-meta.rds")
saveRDS(expr_count, "qc/matrix.rds")

cat(green("successfully saved expression matrix."), crlf)

shared[["counts"]] <- expr_count
shared[["meta_sample_raw"]] <- sample_meta
shared[["meta_gene_raw"]] <- genes_meta
