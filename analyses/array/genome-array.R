
# affymetrix ----------------------------------------------------------------

# 大小写敏感
series_id <- "GSE41627"
download_table <- geo_get_supplemental_links(series_id)
geo_get_supplemental_files(series_id)

expr_set <- geo_get_series_annot(series_id)
platform <- geo_get_platform_annot(expr_set)
groupings <- geo_get_groupings(expr_set)

# 构建文件名称列表
file_names <- paste(groupings[["geo_accession"]], "_",
                    groupings[["title"]], ".cel.gz", sep = "")
groupings <- cbind(groupings, file = file_names)

# 从某个目录下读取一个文件列表的 .cel.gz 文件，返回一个 AffyBatch 对象
read_affy_cel <- function(base_path, targets, filename_col = "file",
                          filename_vec = targets[[filename_col]]) {

  affy::ReadAffy(filenames = filename_vec,
                 celfile.path = base_path,
                 phenoData = targets)
}

affy_batch <- read_affy_cel(base_path = paste("geo", series_id, sep = "/"),
                            targets = groupings,
                            filename_col = "file")
