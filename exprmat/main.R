#!/usr/bin/env Rscript

# construct the shared variable table.

shared <- list()
crlf <- "\n"

# the shared status variable of pipeline.

shared[["is_reference_assigned"]] <- FALSE
shared[["is_loaded"]] <- FALSE
shared[["is_qc"]] <- FALSE
shared[["is_norm"]] <- FALSE
shared[["is_integrate"]] <- FALSE

# load the required libraries

.libPaths(c(
  "/home/yang-z/R/bioinfo/4.4",
  "/usr/lib64/R/library",
  "/usr/share/R/library"
))

suppressPackageStartupMessages({
  require(crayon)
  require(argparse)
  require(dplyr)
  require(stringr)
  require(tibble)
  require(data.table)
  require(purrr)
  require(ggplot2)
  require(extrafont)

  require(SingleCellExperiment)
  require(Seurat)
  require(SeuratWrappers)
  require(scater)
  require(scuttle)

  source("utils.R")
})

# global path resources

gp_scrublet <- "~/Documents/bioinfo/scrublet"
gp_singler <- "~/Documents/bioinfo/singler"
gp_refseq <- "~/Documents/bioinfo/refseq"
gp_annot <- "~/Documents/bioinfo/annot"
gp_base <- "~/Documents/bioinfo/exprmat"

# set up the common repl interface

invoke_command <- function(src, vargs) {
  tryCatch(
    source(paste(gp_base, src, sep = "/")),
    error = function(e) { message(e) }, # nolint
    finally = {}
  )
}

read <- function() {
  cat(crayon::green("$ "))
  readLines("stdin", n = 1, encoding = "utf-8")
}

while (TRUE) { # nolint

  cat(crlf)
  command <- read()
  cat(crlf)

  if (command == "ls") print(names(shared))
  else if (command == "gc") print(gc())
  else if (command == "trace") traceback()
  else if (command == "q") q()
  else if (command == "wd") print(getwd())
  else if (command == ":") cat("call system command with <:>", crlf)
  else if (stringr::str_starts(command, ":")) system(str_sub(command, 2))

  # exprmat-defined function utilities.

  else { # nolint

    vargs <- str_split(command, " ")[[1]]
    if (length(vargs) < 1) {
      next
    } else {

      # extract command target
      cmdtarget <- vargs[1]

      if (cmdtarget %in% c(

        # file system operations
        "cd",

        # reference construction and selection operations
        "gffindex", "gffselect", "gffexonlen", "refer",

        # input source specification
        "read", "read10x", "integrate", "group",

        # quality control
        "qc",

        # normalization
        "norm",

        # advanced single-cell dataset manipulation
        "dimreduc", "cluster", "transfer", "annot", "run",

        # programming language
        "view", "clear", "table",

        # manipulating seurat object
        "chassay", "dim", "intgmeta", "cname"

      )) {

        # extract command arguments
        if (length(vargs) == 1) vargs <- c()
        else vargs <- vargs[2:length(vargs)]
        vargs <- vargs[!is.na(vargs)]
        shared[["vargs"]] <- vargs

        # find and invoke the specified script with the given arguments.
        # arguments will be stored in the shared object.
        # by safety, the sub-scripts should refer only to variables stored in
        # the shared object table. (by convention, except vargs)

        suspected_file_name <- paste(cmdtarget, "R", sep = ".")
        suspected_fpath <- paste(gp_base, suspected_file_name, sep = "/")
        if (file.exists(suspected_fpath))
          invoke_command(suspected_file_name, vargs)

        rm(suspected_file_name, suspected_fpath)
      }
    }

  }

  # here, we will update the shared status after running every command.

  if (file.exists("integrated.rds")) {
    shared[["is_norm"]] <- TRUE
    shared[["is_integrate"]] <- TRUE
  }

  if (file.exists("genome.rds")) {
    shared[["is_reference_assigned"]] <- TRUE
  }

  if (file.exists("features/matrix.rds") &&
        file.exists("features/genes-meta.rds") &&
        file.exists("features/samples-meta.rds")) {
    shared[["is_loaded"]] <- TRUE
  }

  if (file.exists("qc/matrix.rds") &&
        file.exists("qc/genes-meta.rds") &&
        file.exists("qc/samples-meta.rds")) {
    shared[["is_qc"]] <- TRUE
  }

  if (file.exists("norm/linear.rds") &&
        file.exists("norm/log.rds") &&
        file.exists("norm/seurat.rds")) {
    shared[["is_norm"]] <- TRUE
  }
}
