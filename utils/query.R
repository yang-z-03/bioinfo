
geo.links <- function(index) {
  url <- lapply(index, function(x) {
    GEOquery::getGEOSuppFiles(x, makeDirectory = FALSE,
                              fetch_files = FALSE)
  }) |> dplyr::bind_rows()
  return(url)
}

geo.links.wget <- function(index) {
  url <- geo.links(index)
  commands <- paste('wget -O ', url $ fname, ' \'', url $ url, '\'', sep = '')
  cat(commands, sep = '\n')
}
