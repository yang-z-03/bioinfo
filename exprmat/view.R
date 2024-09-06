
parser_view <- argparse::ArgumentParser(prog = "view")

parser_view $ add_argument(type = "character", dest = "object",
                           help = "variable object name to view")

parser_view $ add_argument(
  "-p", "--print", action = "store_const",
  help = "use print method to display object [default]",
  dest = "method", const = "print", default = "print"
)

parser_view $ add_argument(
  "-v", "--view", action = "store_const",
  help = "use view method to display object",
  dest = "method", const = "view"
)

parser_view $ add_argument(
  "-s", "--str", action = "store_const",
  help = "use str method to display object",
  dest = "method", const = "str"
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser_view $ print_help()
  stop()
} else {
  pargs <- parser_view $ parse_args(vargs)
}

if (pargs $ method == "print") {
  print(shared[[pargs $ object]])
} else if (pargs $ method == "view") {
  utils::View(shared[[pargs $ object]])
} else if (pargs $ method == "str") {
  str(shared[[pargs $ object]])
}
