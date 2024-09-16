
# 3.  PSEUDO-TIME AND CELL TRAJECTORIES ========================================

# during development, in response to stimuli, and througout life, cells
# transition from one functional "state" to another. cells in different states
# express different sets of genes, producing a dynamic repetoire of proteins
# and metabolites that carry out their work. as cells move between states,
# undergo a process of transcriptional re-configuration, with some genes being
# silenced and others newly activated. these transient states are often hard to
# characterize because purifying cells in between more stable endpoint states
# can be difficult or impossible. single-cell rna-seq can enable you to see
# these states without the need for purification. however, to do so, we must
# determine where each cell is the range of possible states.


## 3.1  the ordering workflow --------------------------------------------------

# before we get into the details of ordering cells along a trajectory, it's
# important to understand what monocle is doing. the ordering workflow has
# three main steps, each of which involve a significant machine learning task.

# step 1: choosing genes that define progress
#
#   inferring a single-cell trajectory is a machine learning problem. the first
#   step is to select the genes monocle will use as input for its machine
#   learning approach. this is called feature selection, and it has a major
#   impact in the shape of the trajectory. in single-cell rna-seq, genes
#   expressed at low levels are often very noisy, but some contain important
#   information regarding the state of the cell. monocle orders cells by
#   examining the pattern of expression of these genes across the cell
#   population. monocle looks for genes that vary in "interesting"
#   (i.e. not just noisy) ways, and uses these to structure the data. monocle
#   provides you with a variety of tools to select genes that will yield a
#   robust, accurate, and biologically meaningful trajectory. you can use these
#   tools to either perform a completely "unsupervised" analysis, in which
#   monocle has no forehand knowledge of which gene you consider important.
#   alternatively, you can make use of expert knowledge in the form of genes
#   that are already known to define biolgical progress to shape monocle's
#   trajectory. we consider this mode "semi-supervised", because monocle will
#   augment the markers you provide with other, related genes.

# step 2: reducing the dimensionality of the data
#
#   once we have selected the genes we will use to order the cells, monocle
#   applies a dimensionality reduction to the data. monocle uses a recently
#   developed algoithm called reversed graph embedding to reduce the data's
#   dimensionality.

# step 3: ordering the cells in pseudotime
#
#   with the expression data projected into a lower dimensional space, monocle
#   is ready to learn the trajectory that describes how cells transition
#   from one state into another. monocle assumes that the trajectory has a tree
#   structure, with one end of it the "root", and the others the "leaves".
#   monocle's job is to fit the best tree it can to the data. this task is
#   called manifold learning. a cell at the beginning of the biological process
#   starts at the root and progresses along the trunk until it reaches the
#   first branch, if there is one. that cell must then choose a path, and moves
#   further and further along the tree until it reaches a leaf. a cell's
#   pseudotime value is the distance it would have to travel to get back to
#   the root.


## 3.2  selecting genes of trajectories ----------------------------------------

# First, we must decide which genes we will use to define a cell's progress
# through myogenesis. We ultimately want a set of genes that increase (or
# decrease) in expression as a function of progress through the process we're
# studying.

# Ideally, we'd like to use as little prior knowledge of the biology of the
# system under study as possible. We'd like to discover the important ordering
# genes from the data, rather than relying on literature and textbooks, because
# that might introduce bias in the ordering. We'll start with one of the simpler
# ways to do this, but we generally recommend a somewhat more sophisticated
# approach called "dpFeature".


### 3.2.1  temporal experimental designs ---------------------------------------

# you can use any other de packages (e.g. seurat) as you like. it is just
# a designed experiment grouping and de between them.

diff_test_res <- differentialGeneTest(
  expr_myo[expressed_genes, ],
  fullModelFormulaStr = "~ State"
)

ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))


### 3.2.2  select genes based on the difference between clusters ---------------

# we first select superset of feature genes as genes expressed in at least 5%
# of all the cells. (you may achieve this using any other methods)

expr_myo <- detectGenes(expr_myo, min_expr = 0.1)
fData(expr_myo) $ use_for_ordering <-
  fData(expr_myo) $ num_cells_expressed > 0.05 * ncol(expr_myo)

# then we will perform a pca analysis to identify the variance explained by each
# pc (principal component). we can look at a scree plot and determine how many
# pca dimensions are wanted based on whether or not there is a significant gap
# between that component and the component after it.

plot_pc_variance_explained(expr_myo, return_all = FALSE)

# reduce the dimension

expr_myo <- reduceDimension(
  expr_myo, max_components = 2,
  norm_method = "log",
  num_dim = 3, # how many pca dimensions you would like to use?
  reduction_method = "tSNE",
  verbose = TRUE
)

# Then we can run density peak clustering to identify the clusters on the 2-D
# t-SNE space. The densityPeak algorithm clusters cells based on each cell's
# local density (Ρ) and the nearest distance (Δ) of a cell to another cell with
# higher distance. We can set a threshold for the Ρ, Δ and define any cell with
# a higher local density and distance than the thresholds as the density peaks.

# you may run with a default setting first.
expr_myo <- clusterCells(expr_myo, verbose = FALSE)
plot_cell_clusters(expr_myo, color_by = "as.factor(Cluster)")

# and plot the distribution of P and Δ to confirm the threshold.
plot_rho_delta(expr_myo, rho_threshold = 2, delta_threshold = 4)

# and manually give your threshold for P and Δ.
expr_myo <- clusterCells(
  expr_myo,
  rho_threshold = 2,
  delta_threshold = 4,
  skip_rho_sigma = TRUE, # remember to set this. to avoid automatic calc.
  verbose = FALSE
)

# then, we will get the de genes between clusters.

cluster_degenes <- differentialGeneTest(
  expr_myo[expressed_genes, ],
  fullModelFormulaStr = "~ Cluster",
  cores = 1
)

ordering_genes <- row.names(
  cluster_degenes
)[order(cluster_degenes $ qval)][1:1000]


### 3.2.3  selecting highly-variable genes -------------------------------------

# genes that vary a lot are often highly informative for identifying cell
# subpopulations or ordering cells along a trajectory. in rna-seq, a gene's
# variance typically depends on its mean, so we have to be a bit careful about
# how we select genes based on their variance.

disp_table <- dispersionTable(expr_myo)
ordering_genes <- subset(
  disp_table,
  mean_expression >= 0.5 &
    dispersion_empirical >= 1 * dispersion_fit
) $ gene_id


### 3.2.4  supervised (manually-given) markers ---------------------------------

# Unsupervised ordering is desirable because it avoids introducing bias into
# the analysis. However, unsupervised machine learning will sometimes fix on a
# strong feature of the data that's not the focus of your experiment. For
# example, where each cell is in the cell cycle has a major impact on the shape
# of the trajectory when you use unsupervised learning. But what if you wish to
# focus on cycle-independent effects in your biological process? Monocle's
# "semi-supervised" ordering mode can help you focus on the aspects of the
# process you're interested in.

# Ordering your cells in a semi-supervised manner is very simple. You first
# define genes that mark progress using the CellTypeHierchy system, very
# similar to how we used it for cell type classification. Then, you use it to
# select ordering genes that co-vary with these markers. Finally, you order
# the cell based on these genes just as we do in unsupervised ordering.
# So the only difference between unsupervised and semi-supervised ordering is
# in which genes we use for ordering.

# As we saw before, myoblasts begin differnentation by exiting the cell cycle
# and then proceed through a sequence of regulatory events that leads to
# expression of some key muscle-specific proteins needed for contraction. We
# can mark cycling cells with cyclin B2 (CCNB2) and recognize myotubes as
# those cells expressed high levels of myosin heavy chain 3 (MYH3).

ccnb2 <- row.names(subset(fData(expr_myo), gene_short_name == "CCNB2"))
myh3 <- row.names(subset(fData(expr_myo), gene_short_name == "MYH3"))

cth <- newCellTypeHierarchy()

cth <- addCellType(
  cth, "Cycling myoblast",
  classify_func = function(x) {
    x[ccnb2, ] >= 1
  }
)

cth <- addCellType(
  cth, "Myotube",
  classify_func = function(x) {
    x[myh3, ] >= 1
  }
)

cth <- addCellType(
  cth, "Reserve cell",
  classify_func = function(x) {
    x[myh3, ] == 0 & x[ccnb2, ] == 0
  }
)

expr_myo <- classifyCells(expr_myo, cth)

marker_diff <- markerDiffTable(
  expr_myo[expressed_genes, ], cth, cores = 1
)

ordering_genes <-
  row.names(marker_diff)[order(marker_diff $ qval)][1:1000]


## 3.3  later steps ------------------------------------------------------------

# step 1 using selected ordering genes.

expr_myo <- setOrderingFilter(expr_myo, ordering_genes)
plot_ordering_genes(expr_myo)

# next, we will reduce the space down to one with two dimensions, which we will
# be able to easily visualize and interpret while monocle is ordering the cells.

expr_myo <- reduceDimension(expr_myo, max_components = 2, method = "DDRTree")

# finally, calculate the trajectory.

expr_myo <- orderCells(expr_myo)


## 3.4  visualization of trees -------------------------------------------------

# Once the cells are ordered, we can visualize the trajectory in the
# reduced dimensional space.

plot_cell_trajectory(expr_myo, color_by = "Hours")

# The trajectory has a tree-like structure. We can see that the cells collected
# at time zero are located near one of the tips of the tree, while the others
# are distributed amongst the two "branches". Monocle doesn't know a priori
# which of the trajectory of the tree to call the "beginning", so we often
# have to call orderCells again using the root_state argument to specify the
# beginning. First, we plot the trajectory, this time coloring the cells by
# "State":

plot_cell_trajectory(expr_myo, color_by = "State")

# and set the cluster that contains most Hours = 0 cells as the beginning:

state <- function(cds) {
  if (length(unique(pData(cds) $ State)) > 1) {
    t0_counts <- table(pData(cds) $ State, pData(cds) $ Hours)[, "0"]
    return(as.numeric(names(t0_counts)[
      which(t0_counts == max(t0_counts))
    ]))
  } else return(1)
}

expr_myo <- orderCells(expr_myo, root_state = state(expr_myo))
plot_cell_trajectory(expr_myo, color_by = "Pseudotime")

# you can also split the states by separate facets

plot_cell_trajectory(expr_myo, color_by = "State") +
  facet_wrap(~ State, nrow = 1)

# visualize the gene expressions by each state group

blast_genes <- row.names(
  subset(fData(expr_myo),
         gene_short_name %in% c("CCNB2", "MYOD1", "MYOG"))
)

plot_genes_jitter(expr_myo[blast_genes, ], grouping = "State", min_expr = 0.1)

# or by pseudo-time

expressed_genes <-  row.names(
  subset(
    fData(expr_myo),
    num_cells_expressed >= 10
  )
)

filtered <- expr_myo[expressed_genes, ]

my_genes <- row.names(
  subset(
    fData(HSMM_filtered),
    gene_short_name %in% c("CDK1", "MEF2C", "MYH3")
  )
)

cds_subset <- filtered[my_genes, ]
plot_genes_in_pseudotime(cds_subset, color_by = "Hours")
