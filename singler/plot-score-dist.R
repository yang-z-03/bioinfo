
plot_score_dist <- function(
    results,
    show = NULL,
    labels.use = colnames(results $ scores),
    references = NULL,
    scores.use = NULL,
    calls.use = 0,
    pruned.use = NULL,
    size = 0.5,
    ncol = 5,
    dots.on.top = TRUE,
    this.color = "#F0E442",
    pruned.color = "#E69F00",
    other.color = "gray60",
    show.nmads = 3,
    show.min.diff = NULL,
    grid.vars = list()) {

  if (!is.null(show)) {
    show <- match.arg(show, c("scores", "delta.med", "delta.next"))
    if (show != "scores") {
      .Deprecated(new = "plotDeltaDistrbiution")
      return(plot_delta_dist(results,
        show = show, labels.use = labels.use,
        references = scores.use, size = size, ncol = ncol,
        dots.on.top = dots.on.top, this.color = this.color,
        pruned.color = pruned.color, grid.vars = grid.vars
      ))
    } else {
      .Deprecated(old = "show=\"scores\"")
    }
  }

  results <- .ensure_named(results)
  is.combined <- !is.null(results $ orig.results)
  ref.names <- colnames(results $ orig.results)

  if (!is.null(scores.use)) {
    references <- scores.use
    .Deprecated(old = "scores.use", new = "references")
  }
  if (is.null(references)) {
    references <- c(0L, seq_along(results $ orig.results))
  }

  plots <- vector("list", length(references))
  for (i in seq_along(plots)) {

    # Pulling out the scores to use in this iteration.
    chosen <- references[i]
    if (chosen == 0L) {
      current.results <- results
    } else {
      current.results <- results $ orig.results[[chosen]]
    }

    scores <- current.results $ scores
    scores.title <- .values_title(is.combined, chosen, ref.names, show)

    # Pulling out the labels to use in this iteration.
    labels <- current.results $ labels
    labels.title <- .values_title(is.combined, chosen, ref.names, "Labels")

    # Pulling out the pruning calls to use in this iteration.
    prune.calls <- NULL
    if (!is.na(pruned.color)) {
      prune.calls <- current.results $ pruned.labels
    }

    # Actually creating the plot
    plots[[i]] <- .plot_score_distribution(
      scores = scores, labels = labels, prune.calls = prune.calls,
      labels.use = labels.use,
      labels.title = labels.title, scores.title = scores.title,
      this.color = this.color, pruned.color = pruned.color,
      other.color = other.color,
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

.plot_score_distribution <- function(
    scores, labels, prune.calls, labels.use,
    labels.title, scores.title,
    this.color, pruned.color, other.color, size, ncol, dots.on.top) {
  # Create a dataframe with separate rows for each score in values.
  df <- data.frame(
    label = rep(colnames(scores), nrow(scores)),
    values = as.numeric(t(scores)),
    stringsAsFactors = FALSE
  )

  # Add whether this label is the final label given to each cell.
  df $ cell.calls <- rep("other", nrow(df)) # rep() protects when nrow(df)=0.
  is.called <- df $ label == rep(labels, each = ncol(scores))
  df $ cell.calls[is.called] <- "assigned"

  # Replace cell.call if cell was pruned.
  if (!is.null(prune.calls)) {
    is.pruned <- rep(is.na(prune.calls), each = ncol(scores))
    df $ cell.calls[is.pruned & is.called] <- "pruned"
  }

  # Trim dataframe by labels
  keep <- df $ label %in% labels.use
  if (any(keep)) {
    df <- df[keep, ]
  } else {
    warning("ignoring 'labels.use' as it has no values in ", scores.title)
  }

  # Making the violin plots.
  p <- ggplot2::ggplot(
    data = df,
    ggplot2::aes_string(x = "cell.calls", y = "values", fill = "cell.calls")
  ) +
    ggplot2::scale_fill_manual(
      name = labels.title,
      breaks = c("assigned", "pruned", "other"),
      values = c(this.color, pruned.color, other.color)
    )

  jit <- ggplot2::geom_jitter(
    height = 0, width = 0.3, color = "black",
    shape = 16, size = size, na.rm = TRUE
  )

  .pretty_violins(p,
    df = df, ncol = ncol, scores.title = scores.title,
    size = size, dots.on.top = dots.on.top, jitter = jit
  )
}

.pretty_violins <- function(
    p, df, ncol, scores.title, size, dots.on.top, jitter, ...
  ) {

  p <- p + ggplot2::theme_classic() +
    ggplot2::facet_wrap(facets = ~label, ncol = ncol) +
    ggplot2::ylab(scores.title)

  if (nlevels(as.factor(df $ label)) == 1) {
    p <- p + ggplot2::scale_x_discrete(name = NULL, labels = NULL)
  } else {
    p <- p + ggplot2::scale_x_discrete(name = "Labels", labels = NULL)
  }

  if (!dots.on.top) {
    p <- p + jitter
  }

  p <- p + ggplot2::geom_violin(na.rm = TRUE, ...)

  if (dots.on.top) {
    p <- p + jitter
  }

  p
}
