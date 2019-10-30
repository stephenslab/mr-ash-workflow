
#'
#'
#'
#' Plotting routines

## function for boxplot
my.box <- function (dat, x, y, cols = gg_color_hue(13),
                    shapes = 1:15) {
  return(ggplot(dat,aes_string(x = x, y = y, fill = "fit")) +
           scale_shape_manual(values = shapes) +
           geom_boxplot(alpha = 0.1, aes(color = fit), outlier.alpha = 0) +
           scale_alpha_manual(values = 0.1) +
           scale_fill_manual(values = cols, guide = "none") +
           scale_color_manual(values = cols, guide = "none") +
           labs(x           = "") +
           theme_cowplot(font_size = 14))
}

## function for boxplot2
my.box2 <- function (dat, x, y, cols = gg_color_hue(13),
                    shapes = 1:15) {
  return(ggplot(dat,aes_string(x = x, y = y, fill = "fit")) +
           scale_shape_manual(values = shapes) +
           stat_summary(fun.y= mean, colour= cols, geom="point", 
                        shape = 18, size = 3) + 
           geom_boxplot(alpha = 0.1, aes(color = fit), outlier.alpha = 0) +
           scale_alpha_manual(values = 0.1) +
           scale_fill_manual(values = cols, guide = "none") +
           scale_color_manual(values = cols, guide = "none") +
           labs(x           = "") +
           theme_cowplot(font_size = 14))
}


## function for line
my.line <- function (dat, x, y, cols = NULL,
                     shapes = NULL) {
  return(ggplot(dat,aes_string(x = x, y = y, fill = "fit")) +
           scale_shape_manual(values = shapes) +
           scale_color_manual(values = cols) +
           geom_point(aes(color = fit, shape = fit)) +
           geom_line(aes(color = fit)) +
           scale_fill_discrete(guide = "none") +
           labs(x           = "") +
           theme_cowplot(font_size = 14))
}

## function for violin plot
my.jitter <- function (dat, x, y,
                       values = c(1,2,0,3,4,5)) {
  return(ggplot(dat,aes_string(x = x, y = y, fill = "fit")) +
           geom_jitter(aes(color = fit, shape = fit)) +
           scale_shape_manual(values = values) +
           geom_violin(alpha = 0.1, aes(color = fit, fill = fit), scale = "width") +
           scale_alpha_manual(0.1) +
           scale_fill_discrete(guide = "none") +
           labs(x           = "") +
           theme_cowplot(font_size = 14))
}

## function for color
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
