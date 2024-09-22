
parser <- argparse::ArgumentParser(prog = "view")

parser $ add_argument(type = "character", dest = "object",
                      help = "variable object name to view")

parser $ add_argument(
  "-p", "--print", action = "store_const",
  help = "use print method to display object [default]",
  dest = "method", const = "print", default = "print"
)

parser $ add_argument(
  "-v", "--view", action = "store_const",
  help = "use view method to display object",
  dest = "method", const = "view"
)

parser $ add_argument(
  "-s", "--str", action = "store_const",
  help = "use str method to display object",
  dest = "method", const = "str"
)

parser $ add_argument(
  "-c", dest = "cols", type = "character",
  default = c(), nargs = "+",
  help = "see only specified column (only work for tibbles)"
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

if (pargs $ method == "print") {
  print(shared[[pargs $ object]])
} else if (pargs $ method == "view") {
  utils::View(shared[[pargs $ object]])
} else if (pargs $ method == "str") {
  str(shared[[pargs $ object]])
}

cat(crlf)
objclass <- shared[[pargs $ object]] |> class()
cat(red("class:"), objclass, crlf)

switch(
  objclass[1],
  Seurat = {
    cat(yellow("metadata columns:"), crlf)
    print(shared[[pargs $ object]] @ meta.data |> colnames())
    cat(yellow("gene names"), crlf)
    print(shared[[pargs $ object]] |> rownames() |> head(20))
    cat(yellow("sample names"), crlf)
    print(shared[[pargs $ object]] |> colnames() |> head(20))
  },
  tbl_df = {
    cat(yellow("full column names"), crlf)
    if (shared[[pargs $ object]] |> colnames() |> length() <= 200) {
      print(shared[[pargs $ object]] |> colnames())
    } else {
      cat(yellow("  more than 200 columns, omitted."), crlf)
    }

    if (pargs $ cols |> length() > 0) {
      cat(yellow("selected columns view:"), crlf)
      print(shared[[pargs $ object]][, pargs $ cols])
    }
  }
)