
if (!shared[["is_norm"]]) {
  cat(red("you should generate a dataset before using analyses function 'w'."))
  cat(crlf)
  stop()
}

cat(blue("saving gene metadata ..."), crlf)
saveRDS(shared[["meta_gene"]], "norm/genes-meta.rds")
cat(blue("saving sample metadata ..."), crlf)
saveRDS(shared[["meta_sample"]], "norm/samples-meta.rds")
cat(blue("saving seurat object ..."), crlf)
saveRDS(shared[["seurat"]], "norm/seurat.rds")
