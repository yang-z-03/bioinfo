#!/usr/bin/env Rscript

# construct the shared variable table.

shared <- list()
crlf <- "\n"

# the shared status variable of pipeline.

shared[["is_reference_assigned"]] <- FALSE
shared[["is_loaded"]] <- FALSE
shared[["is_qc"]] <- FALSE
shared[["is_norm"]] <- FALSE
shared[["is_ready"]] <- FALSE

parser <- argparse::ArgumentParser(
  prog = "exprmat",
  description = "pipeline for expression matrix"
)

parser $ add_argument(
  "-s", type = "character", dest = "script", default = "",
  help = "the input script file, will be executed line by line"
)

parser $ add_argument(
  "--local-libs", dest = "local.lib", default = FALSE, action = "store_true",
  help = "use local library"
)

pargs <- parser $ parse_args()

# here, we resolve the line cutters

parse_script <- function(fname) {
  if (fname |> stringr::str_trim() |> stringr::str_length() == 0)
    return(c())

  cmdline <- NULL
  if (file.exists(fname)) {
    cmdline <- readLines(con = file(fname, "r"), n = -1)
  } else {
    cat(crayon::red(paste(fname, "not found.")))
    return(c())
  }

  # here, we resolve the line cutters
  cmdlist <- NULL
  appendline <- FALSE
  for (line in cmdline) {
    if (appendline) {
      cmdlist[length(cmdlist)] <- paste(
        cmdlist[length(cmdlist)] |>
          substr(1, stringr::str_length(cmdlist[length(cmdlist)]) - 1),
        line |> stringr::str_trim(), sep = " "
      )
    } else {
      cmdlist <- c(cmdlist, line |> stringr::str_trim())
    }

    if (line |> stringr::str_trim() |> stringr::str_ends("\\\\")) {
      appendline <- TRUE
    } else {
      appendline <- FALSE
    }
  }

  return(cmdlist)
}

cmdlist <- parse_script(pargs $ script)

if (pargs $ local.lib) {
  .libPaths(c(
    "/home/yang-z/R/bioinfo/4.4",
    "/usr/lib64/R/library",
    "/usr/share/R/library"
  ))
}

# global path resources

gp_scrublet <- "~/bioinfo/scrublet"
gp_singler <- "~/bioinfo/singler"
gp_refseq <- "~/bioinfo/refseq"
gp_gencode <- "~/bioinfo/gencode"
gp_annot <- "~/bioinfo/annot"
gp_base <- "~/bioinfo/exprmat"

setwd(gp_base)

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
  require(reticulate)
  reticulate::use_python("/home/data/yangzhen/miniconda/bin/python")

  require(SingleCellExperiment)
  require(SC3)
  require(Seurat)
  require(SeuratWrappers)
  require(scater)
  require(scuttle)

  source("utils.R")
})

# set up the common repl interface

invoke_command <- function(src, vargs) {
  tryCatch(
    source(paste(gp_base, src, sep = "/")),
    error = function(e) { message(e) }, # nolint
    finally = {}
  )
}

auto_mode <- FALSE
read <- function() {
  cwd <- basename(getwd())
  if (auto_mode) cat(crayon::red("(auto) "))
  cat(crayon::yellow(cwd))
  cat(crayon::green("$ "))
  readLines("stdin", n = 1, encoding = "utf-8")
}

shlex <- reticulate::import("shlex")

while (TRUE) { # nolint

  cat(crlf)
  if (length(cmdlist) == 0) {
    auto_mode <- FALSE
    command <- read()
  } else {
    auto_mode <- TRUE
    command <- cmdlist[1]
    cmdlist <- cmdlist[-1]

    if (auto_mode) cat(crayon::red("(auto) "))
    cat(crayon::yellow(basename(getwd())))
    cat(crayon::green("$ "))
    cat(command)

    # detect comment chars
    cmt <- str_locate(command, "#")[1, "start"] - 1
    if (!is.na(cmt)) command <- substr(command, 1, cmt)

    cat(crlf)
  }
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

    # shell split

    vargs <- shlex $ split(command)

    if (length(vargs) < 1) {
      next
    } else {

      # extract command target
      cmdtarget <- vargs[1]

      if (cmdtarget %in% c(

        # file system operations
        "cd",

        # reference construction and selection operations
        "gffindex", "gffselect", "gffexonlen", "refer", "refer.gencode",

        # input source specification
        "read", "read10x", "integrate", "group", "intgmeta", "intsmeta",

        # quality control
        "qc",

        # normalization
        "norm",

        # advanced single-cell dataset manipulation
        "dimreduc", "cluster", "annot", "run", "de", "lsclust", "rnclust",
        "cpsclust", "defclust", "addclust", "cpdbin",

        # programming language
        "view", "clear", "table",

        # manipulating seurat object
        "chassay", "dim", "cname", "w",

        # proteomics
        "readp", "normp",

        "source"

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

  if (file.exists("norm/genes-meta.rds") &&
      file.exists("norm/samples-meta.rds") &&
      file.exists("norm/seurat.rds")) {
    shared[["is_ready"]] <- TRUE
  }
}
