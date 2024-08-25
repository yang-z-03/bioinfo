
library(fs)

root_dir <- paste("geo", "GSE139324", sep = "/")

arrange_files_to_10x <- function(root_dir) {
  files <- list.files(root_dir)
  matrix <- strsplit(files, split = "_")

  for (x in matrix) {
    if (length(x) != 5) next

    attempt_dir0 <- paste(root_dir, paste(x[1], x[2], x[3], sep = "_"), sep = "/")
    attempt_dir <- paste(root_dir, paste(x[1], x[2], x[3], sep = "_"), x[4], sep = "/")
    if (!dir.exists(attempt_dir0)) {
      dir.create(attempt_dir0)
    }
    if (!dir.exists(attempt_dir)) {
      dir.create(attempt_dir)
    }

    file_move(
      paste(root_dir, paste(x[1], x[2], x[3], x[4], x[5], sep = "_"), sep = "/"),
      paste(root_dir, paste(x[1], x[2], x[3], sep = "_"), x[4], x[5], sep = "/")
    )
  }
}

arrange_files_to_10x(root_dir)
