
plot_delta_dist <- function(
    results,
    show = c("delta.med", "delta.next"),
    labels.use = colnames(results $ scores),
    references = NULL,
    chosen.only = TRUE,
    size = 2,
    ncol = 5,
    dots.on.top = TRUE,
    this.color = "#000000",
    pruned.color = "#E69F00",
    grid.vars = list()
  ) {

  results <- .ensure_named(results)
  show <- match.arg(show)
  is.combined <- !is.null(results $ orig.results)
  ref.names <- colnames(results $ orig.results)

  if (is.null(references)) {
    if (is.combined) {
      # Combined 'delta.med'/'delta.next' are not worth showing.
      references <- seq_along(results $ orig.results)
    } else {
      references <- 0L
    }
  }

  plots <- vector("list", length(references))
  for (i in seq_along(plots)) {
    # Pulling out the scores to use in this iteration.
    chosen <- references[i]
    if (is.combined) {
      if (chosen == 0L) {
        stop("deltas cannot be shown for combined results")
      }
      current.results <- results $ orig.results[[chosen]]
      if (chosen.only) {
        current.results <- current.results[chosen == results $ reference, ]
      }
    } else {
      current.results <- results
    }
    scores.title <- .values_title(is.combined, chosen, ref.names, show)

    # Computing the values that we want to show.
    if (show == "delta.med") {
      values <- delta_from_median(current.results)
    } else if (show == "delta.next") {
      values <- current.results $ delta.next
    }

    # Pulling out the labels to use in this iteration.
    labels <- current.results $ labels
    labels.title <- .values_title(is.combined, chosen, ref.names, "Labels")

    # Checking if we need pruning.
    pruned <- NULL
    if (!is.na(pruned.color)) {
      pruned <- is.na(current.results $ pruned.labels)
    }

    # Actually creating the plot
    plots[[i]] <- .plot_delta_distribution(
      values = values, labels = labels, pruned = pruned, 
      labels.use = labels.use,
      labels.title = labels.title, scores.title = scores.title,
      this.color = this.color, pruned.color = pruned.color,
      size = size, ncol = ncol, dots.on.top = dots.on.top
    )
  }

  if (length(plots) == 1L) {
    # Doing this to be consistent with raw ggplot output.
    plots[[1]]
  } else {
    if (!is.null(grid.vars) && length(references) > 1L) {
      do.call(gridExtra::grid.arrange, c(plots, grid.vars))
    } else {
      plots
    }
  }
}

.plot_delta_distribution <- function(
    values, labels, pruned, labels.use,
    labels.title, scores.title,
    this.color, pruned.color,
    size, ncol, dots.on.top) {
  df <- data.frame(
    values = values,
    label = labels,
    x = character(length(values))
  )

  aes.jit <- NULL
  if (!is.null(pruned)) {
    df $ pruned <- pruned
    aes.jit <- ggplot2::aes_string(color = "pruned")
  }

  # Trim dataframe by labels:
  if (any(keep <- labels.use %in% df $ label)) {
    labels.use <- labels.use[keep]
    df <- df[df $ label %in% labels.use, ]
  } else {
    warning("ignoring 'labels.use' as it has no values in ", scores.title)
  }

  # Making the violin plots.
  p <- ggplot2::ggplot(data = df, ggplot2::aes_string(x = "x", y = "values")) +
    ggplot2::xlab("")

  if (!is.null(pruned)) {
    p <- p + ggplot2::scale_color_manual(
      name = "Pruned",
      breaks = c("FALSE", "TRUE"),
      values = c(this.color, pruned.color)
    )
  }

  jit <- ggplot2::geom_jitter(
    mapping = aes.jit, height = 0, width = 0.3,
    shape = 16, size = size, na.rm = TRUE
  )

  .pretty_violins(p,
    df = df, ncol = ncol, scores.title = scores.title,
    size = size, dots.on.top = dots.on.top, jitter = jit, fill = "grey"
  )
}
