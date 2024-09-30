
parser <- argparse::ArgumentParser(prog = "w")

parser $ add_argument("-o", type = "character", dest = "object", default = "",
                      help = "write object, press . for seurat objects")

parser $ add_argument("-d", type = "character", dest = "f", default = "obj.rds",
                      help = "dump file path (better with extension .rds)")

parser $ add_argument("-u", dest = "compress",
                      default = TRUE, action = "store_false",
                      help = "do not compress the rds file")

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

if (pargs $ object != ".") {

  cat(
    blue("saving object of class"),
    red(shared[[pargs $ object]] |> class() |> paste()), crlf
  )
  saveRDS(shared[[pargs $ object]], file = pargs $ f,
          compress = pargs $ compress)

} else {

  if (!shared[["is_ready"]]) {
    cat(red("you should generate a dataset before using analyses function 'w'"))
    cat(crlf)
    stop()
  }

  cat(blue("saving gene metadata ..."), crlf)
  saveRDS(shared[["meta_gene"]], "norm/genes-meta.rds",
          compress = pargs $ compress)
  cat(blue("saving sample metadata ..."), crlf)
  saveRDS(shared[["meta_sample"]], "norm/samples-meta.rds",
          compress = pargs $ compress)
  cat(blue("saving seurat object ..."), crlf)
  saveRDS(shared[["seurat"]], "norm/seurat.rds",
          compress = pargs $ compress)
}