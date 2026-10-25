# Script used to create the package hex sticker (not part of the package).

library(plot3D)
library(hexSticker)

# A Swiss roll as a surface, colored along the roll
t <- seq(1.5 * pi, 4.5 * pi, length.out = 300)
h <- seq(-7, 7, length.out = 40)
grid <- mesh(t, h)

png("data-raw/swiss-roll.png", width = 2000, height = 2000, res = 300, bg = "transparent")
par(mar = c(0, 0, 0, 0))
surf3D(
  x = grid$x * sin(grid$x),
  y = grid$y,
  z = grid$x * cos(grid$x),
  colvar = grid$x,
  col = rainbow(256, end = 0.8),
  border = NA,
  shade = 0.35,
  theta = 20,
  phi = 20,
  box = FALSE,
  colkey = FALSE
)
dev.off()

# Trim the transparent margins
roll <- magick::image_trim(magick::image_read("data-raw/swiss-roll.png"))
magick::image_write(roll, "data-raw/swiss-roll.png")

sysfonts::font_add_google("Roboto Condensed", "roboto")
showtext::showtext_auto()

sticker(
  subplot = "data-raw/swiss-roll.png",
  package = "manifoldsculpting",
  p_family = "roboto",
  p_fontface = "bold",
  p_size = 26,
  p_y = 1.48,
  p_color = "white",
  s_x = 1,
  s_y = 0.86,
  s_width = 0.46,
  h_fill = "black",
  h_color = "#9085e9",
  h_size = 1.4,
  url = "https://github.com/ChristianGoueguel/manifold-sculpting",
  u_color = "white",
  u_size = 5,
  u_family = "roboto",
  spotlight = TRUE,
  l_x = 1,
  l_y = 0.86,
  l_width = 3,
  l_height = 3,
  l_alpha = 0.15,
  filename = "data-raw/sticker.png",
  dpi = 600
)

usethis::use_logo(img = "data-raw/sticker.png", geometry = "240x278", retina = TRUE)

