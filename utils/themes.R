
library(ggplot2)

unify_theme_font <- function(base_size = 11, base_family = "SF Pro Text") {

  txt <- element_text(size = base_size, colour = "black", face = "plain")
  bold_txt <- element_text(size = base_size, colour = "black", face = "bold")

  theme(
    legend.key = element_blank(),
    strip.background = element_blank(),
    text = txt,
    plot.title = txt,
    axis.title = txt,
    axis.text = txt,
    legend.title = bold_txt,
    legend.text = txt
  )
}
