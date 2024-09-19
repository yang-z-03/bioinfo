
# 5.  VISUALIZE TEMPORAL EXPRESSION PATTERNS ===================================

# a common question that arises when studying time-series gene expression
# studies is: "which genes follow similar kinetic trends"? monocle can help
# you answer this question by grouping genes that have similar trends, so you
# can analyze these groups to see what they have in common. monocle provides a
# convenient way to visualize all pseudotime-dependent genes. the function
# plot_pseudotime_heatmap takes a celldataset object (usually containing a only
# subset of significant genes) and generates smooth expression curves much like
# plot_genes_in_pseudotime. then, it clusters these genes and plots them using
# the pheatmap package. this allows you to visualize modules of genes that
# co-vary across pseudotime.

diff_test_res <- differentialGeneTest(
  expr_myo[marker_genes, ],
  fullModelFormulaStr = "~ sm.ns(Pseudotime)"
)

sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
# selected genes that correlates with pseudo time.

plot_pseudotime_heatmap(
  expr_myo[sig_gene_names, ],
  num_clusters = 3, # clustering the genes, by their likely temporal pattern.
  cores = 1,
  show_rownames = TRUE
)