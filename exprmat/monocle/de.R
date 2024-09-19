
# 4.  FINDING GENES THAT CHANGE AS A FUNCTION OF PSEUDO-TIME ===================

# monocle's main job is to put cells in order of progress through a biological
# process (such as cell differentiation) without knowing which genes to look at
# ahead of time. once it's done so, you can analyze the cells to find genes
# that changes as the cells make progress. for example, you can find genes that
# are significantly upregulated as the cells "mature". let's look at a panel of
# genes important for myogenesis:

to_be_tested <- row.names(
  subset(fData(expr),
         gene_short_name %in% c("MYH3", "MEF2C", "CCNB2", "TNNT1"))
)
cds_subset <- expr_myo[to_be_tested, ]

# again, we'll need to specify the model we want to use for differential
# analysis. this model will be a bit more complicated than the one we used to
# look at the differences between celltype. monocle assigns each cell a
# "pseudotime" value, which records its progress through the process in the
# experiment. the model can test against changes as a function of this value.
# monocle uses the vgam package to model a gene's expression level as a smooth,
# nonlinear function of pseudotime.

diff_test_res <- differentialGeneTest(
  cds_subset, fullModelFormulaStr = "~ sm.ns(Pseudotime)"
)

# the sm.ns function states that monocle should fit a natural spline through
# the expression values to help it describe the changes in expression as a
# function of progress. we'll see what this trend looks like in just a moment.
# other smoothing functions are available. once again, let's add in the gene
# annotations so it's easy to see which genes are significant.

diff_test_res[, c("gene_short_name", "pval", "qval")]

# we can plot the expression levels of these genes, all of which show
# significant changes as a function of differentiation, using the function
# plot_genes_in_pseudotime. this function has a number of cosmetic options you
# can use to control the layout and appearance of your plot.

plot_genes_in_pseudotime(cds_subset, color_by = "Hours")
