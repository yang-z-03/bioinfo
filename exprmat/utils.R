
make_matrix <- function(df, rownames = NULL) {
  my_matrix <-  as.matrix(df)
  if (!is.null(rownames))
    rownames(my_matrix) <- rownames
  my_matrix
}

unify_theme_font <- function(base_size = 13, base_family = "Myriad Pro Cond") {
  txt <- ggplot2::element_text(
    colour = "black", face = "plain",
    family = base_family
  )

  bold_txt <- ggplot2::element_text(
    colour = "black", face = "bold",
    family = base_family
  )

  ggplot2::theme_bw(base_size = base_size) + ggplot2::theme(
    legend.key = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    text = txt,
    plot.title = txt,
    axis.title = txt,
    axis.text = txt,
    legend.title = bold_txt,
    legend.text = txt,
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      color = "#e0e0e0",
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(), # remove minor grids
    axis.ticks.length = ggplot2::unit(-0.05, "in") # ticks inwards
  )
}
