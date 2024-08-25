
rebuild_index <- function(trained, num.threads = 1) {
  if (.is_solo(trained)) {
    trained <- .rebuild_index(trained, num.threads = num.threads)
  } else {
    for (i in seq_along(trained)) {
      trained[[i]] <- .rebuild_index(
        trained[[i]],
        num.threads = num.threads
      )
    }
  }
  trained
}

.rebuild_index <- function(trained, num.threads) {
  if (!is(trained$built, "externalptr") || 
      !is_valid_built(trained$built)) {

    trained $ built <- .build_index(
      ref = trained $ ref,
      markers = trained $ markers $ full,
      labels = trained $ labels $ full,
      ulabels = trained $ labels $ unique,
      approximate = trained $ options $ approximate,
      num.threads = num.threads
    )
  }
  trained
}
