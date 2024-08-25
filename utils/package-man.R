
set_repository <- function() {

  # setting the R package mirrors
  options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  options("BioC_mirror" = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor/")

  # display the current mirror site origins of these package sources.
  option <- options()
  option[["repos"]][["CRAN"]]
  option[["BioC_mirror"]]

  # the bioconductor repositories
  BiocManager::repositories()

  # returns the CRAN and Bioconductor mirror sites set.
  return(c(option[["repos"]][["CRAN"]], option[["BioC_mirror"]]))
}

install_packages <- function(pkgs, j = 4) {
  set_repository()

  # first install bioc-manager, skip if already installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", ask = FALSE, update = FALSE)
  }

  # install stringr
  if (!requireNamespace("stringr", quietly = TRUE)) {
    install.packages("stringr", ask = FALSE, update = FALSE)
  }

  # remove the packages that occurs more than once. also detect names that seems
  # to be github addresses.
  pkgs <- unique(pkgs)
  pkgs2 <- pkgs
  logi <- stringr::str_detect(pkgs2, "/")
  pkgs2[logi] <- stringr::str_match(pkgs2[logi], ".*/(.*)$")[, 2]

  # install the packages in the list if currently have not already installed
  new <- !(sapply(pkgs2, requireNamespace, quietly = TRUE))

  if (sum(new) > 0) {
    # display the package newly installed here.
    cat("packages to install: ", pkgs[new], "\n")
  } else {
    cat("all packages given have already installed \n")
  }

  # install packages
  if (any(new)) BiocManager::install(pkgs[new], ask = FALSE, update = FALSE,
                                     Ncpus = j) # num of threads to compile.
}

# the location to the windows compiled binary package, commonly a .zip archive
install_local_binary <- function(name, location) {
  if (!requireNamespace(name, quietly = TRUE)) {
    install.packages(location, repos = NULL, type = "win.binary")
  }
}

# the location to the source code, commonly a .tar.gz archive, may also be a
# .zip archive downloaded from github page.
install_local_source <- function(name, location) {
  if (!requireNamespace(name, quietly = TRUE)) {
    install.packages(location, repos = NULL, type = "source")
  }
}

packs_require <-
  c("tidyverse", "limma", "affy", "oligo", "lumi",
    "beadarray", "GEOquery", "simpleaffy", "gcrma", "readxl",
    "impute", "genefilter", "pd.hugene.1.0.st.v1", "pd.hg.u133.plus.2",
    "AnnotationDbi", "org.Hs.eg.db",
    "hgug4112a.db", "AgiMicroRna", "sva", "DESeq2", "edgeR",
    "lumiHumanIDMapping", "remotes", "pheatmap", "shiny", "aggregation",
    "tidyverse/dplyr", "limma", "hwriter", "devtools", "ggplot2",
    "Seurat", "RColorBrewer", "glmGamPoi", "scales", "SingleCellExperiment",
    "scater", "scran", "patchwork", "reticulate", "data.tree", "RJSONIO")

get_installed_package_list <- function() {
  info <- installed.packages()
  less <- info[, c("Package", "Version", "License", "Built")]
  rownames(less) <- c()
  return(less)
}

packs <- get_installed_package_list()
# install_packages(packs_require)
