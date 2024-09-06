
make_matrix <- function(df, rownames = NULL){
  my_matrix <-  as.matrix(df)
  if (!is.null(rownames))
    rownames(my_matrix) <- rownames
  my_matrix
}