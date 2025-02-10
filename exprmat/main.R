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
  "--conda", type = "character", dest = "conda", default = "base",
  help = "use miniconda environment"
)

parser $ add_argument(
  "--dry-run", dest = "dry.run", default = FALSE, action = "store_true",
  help = "prompt commands but do not run"
)

main_pargs <- parser $ parse_args()

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
    } else if (line |> stringr::str_trim() |> stringr::str_length() > 0) {
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

cmdlist <- parse_script(main_pargs $ script)

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

  require(SingleCellExperiment)
  require(SC3)
  require(Seurat)
  require(SeuratWrappers)
  require(scater)
  require(scuttle)

  source("utils.R")
  reticulate::use_condaenv(main_pargs $ conda)
})

# set up the common repl interface

invoke_command <- function(src, vargs) {
  tryCatch(
    source(paste(gp_base, src, sep = "/")),
    error = function(e) { message(e) }, # nolint
    finally = {}
  )
}

read <- function() {
  cwd <- basename(getwd())
  cat(crayon::yellow(cwd))
  cat(crayon::green(" $ "))
  readLines("stdin", n = 1, encoding = "utf-8")
}

shlex <- reticulate::import("shlex")

.current_index <- 1 # the line of code executing.
.loop_restores <- c() # the line of loop startings
.loop_break <- NULL # the line of loop ending
.loop_indices <- c()
.is_redirected <- FALSE

.for_parser <- argparse::ArgumentParser(
  prog = "for", description = "program control loops"
)

.for_parser $ add_argument(
  "--in", type = "character", dest = "iter",
  help = paste("index-able object, which should return true for either",
               "is.data.frame or is.atomic")
)

.for_parser $ add_argument(
  "ar", type = "character", nargs = "+",
  help = paste("iterating objects")
)

command <- ""

system("clear")

while (TRUE) { # nolint

  if (command != "" &&
        (!main_pargs $ dry.run) &&
        (!.is_control)) cat(crlf)

  .is_control <- FALSE

  # the cmdlist variable stores the input command with automatic execution
  # here, we only implement the most-used controls: the loops.
  # we do not aim to be complete script language ...

  .cmdlen <- length(cmdlist)

  if (.current_index > .cmdlen) {
    # the script is running over. read from the input.
    cat(crayon::green(.current_index), "") # line no.
    command <- read()
    cmdlist <- c(cmdlist, command)
    .current_index <- length(cmdlist)

  } else {
    # read from the cmdlist
    commandstr <- cmdlist[.current_index] |> stringr::str_trim()
    # detect comment chars
    cmt <- str_locate(commandstr, "#")[1, "start"] - 1

    command <- commandstr
    if (!is.na(cmt)) command <- substr(commandstr, 1, cmt)

    # by now, we need to detect for loops. we will parse the sentence.
    vargs <- shlex $ split(command)

    if (vargs[1] == "for") {
      # the loop starter.
      for_stmt <- .for_parser $ parse_args(vargs[-1])
      .is_control <- TRUE

      iterator <- shared[[for_stmt $ iter]]
      validiter <- is.data.frame(iterator) | is.atomic(iterator)
      if (!validiter) {
        cat(crayon::red("fatal: for iterator not indexable!"), crlf)
        stop()
      }

      if (!.is_redirected) { # is not redirected for
        # the initialization of a for loop.
        .loop_restores <- c(.loop_restores, .current_index)
        .loop_indices <- c(.loop_indices, 1)
      } else .is_redirected <- FALSE

      hasitem <- FALSE
      current_loop_index <- .loop_indices |> dplyr::last()
      if (is.data.frame(iterator)) {
        hasitem <- nrow(iterator) >= current_loop_index
      } else if (is.atomic(iterator)) {
        hasitem <- length(iterator) >= current_loop_index
      }

      if (!hasitem) {
        # exit the loop
        .loop_restores <- .loop_restores[1 : (length(.loop_restores) - 1)]
        .loop_indices <- .loop_indices[1 : (length(.loop_indices) - 1)]
        .current_index <- .loop_break
        .loop_break <- NULL

      } else {
        if (is.data.frame(iterator)) {
          coliter <- colnames(iterator)
          for (coli in seq_along(coliter)) {
            shared[[paste("_", for_stmt $ ar[coli], sep = "")]] <-
              iterator[current_loop_index, coliter[coli]]
          }
        } else if (is.atomic(iterator)) {
          shared[[paste("_", for_stmt $ ar[1], sep = "")]] <-
            iterator[current_loop_index]
        }
      }

    } else if (vargs[1] == "endfor") {
      .is_redirected <- TRUE
      .is_control <- TRUE
      .loop_break <- .current_index
      .loop_indices[length(.loop_indices)] <-
        .loop_indices[length(.loop_indices)] + 1
      .current_index <- .loop_restores |> dplyr::last() - 1

    } else if (vargs[1] == "dry") {
      .is_control <- TRUE
    } else if (vargs[1] == "wet") {
      .is_control <- TRUE
    } else if (vargs[1] == "echo") {
      .is_control <- TRUE
      cat(crayon::green(.current_index), "") # line no.
      # now substitute the command strings.
      varnames <- names(shared)
      varnames <- varnames[stringr::str_starts(varnames, "_")]
      varnames <- stringr::str_sub(varnames, 2)
      
      for (v in varnames) {
        command <- gsub(
          paste("\\$\\{", v, "\\}", sep = ""),
          paste("\"", shared[[paste("_", v, sep = "")]], "\"", sep = ""),
          command
        )
        
        command <- gsub(
          paste("\\$\\(", v, "\\)", sep = ""),
          shared[[paste("_", v, sep = "")]],
          command
        )
      }
      
      cat(command |> stringr::str_sub(6)) # remove the echo word.
      cat(crlf)
      
    } else {
      # for running sentense, we should print it out
      cat(crayon::green(.current_index), "") # line no.
      cat(crayon::red("(auto) "))
      cat(crayon::yellow(basename(getwd())))
      cat(crayon::green(" $ "))

      # now substitute the command strings.
      varnames <- names(shared)
      varnames <- varnames[stringr::str_starts(varnames, "_")]
      varnames <- stringr::str_sub(varnames, 2)

      for (v in varnames) {
        commandstr <- gsub(
          paste("\\$\\{", v, "\\}", sep = ""),
          paste("\"", shared[[paste("_", v, sep = "")]], "\"", sep = ""),
          commandstr
        )

        commandstr <- gsub(
          paste("\\$\\(", v, "\\)", sep = ""),
          shared[[paste("_", v, sep = "")]],
          commandstr
        )

        command <- gsub(
          paste("\\$\\{", v, "\\}", sep = ""),
          paste("\"", shared[[paste("_", v, sep = "")]], "\"", sep = ""),
          command
        )

        command <- gsub(
          paste("\\$\\(", v, "\\)", sep = ""),
          shared[[paste("_", v, sep = "")]],
          command
        )
      }

      cat(commandstr)
      cat(crlf)
    }
  }

  # RUN -----------------------------------------------------------------------

  if (main_pargs $ dry.run) {
    .current_index <- .current_index + 1
    if (command == "q") q()
    else if (command == "wet") main_pargs $ dry.run <- FALSE
    next
  }

  if (command != "" && (!.is_control)) cat(crlf)

  if (command == "ls") print(names(shared))
  else if (command == "gc") print(gc())
  else if (command == "q") q()
  else if (command == "wd") print(getwd())
  else if (command == "dry") main_pargs $ dry.run <- TRUE
  else if (command == ":") cat("call system command with <:>", crlf)
  else if (stringr::str_starts(command, ":")) system(str_sub(command, 2))

  # exprmat-defined function utilities.

  else { # nolint

    vargs <- shlex $ split(command) # shell split

    if (length(vargs) < 1) {
      # do nothing
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
        "cpsclust", "defclust", "addclust", "cpdbin", "cpdbhom",

        # programming language
        "view", "clear", "table",

        # manipulating seurat object
        "chassay", "dim", "cname", "w",

        # proteomics
        "readp", "normp",

        "source", "set"

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

  .current_index <- .current_index + 1

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
