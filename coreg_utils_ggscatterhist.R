



my_ggscatterhist <- function (data, x, y, group = NULL, color = "black", fill = NA, 
          palette = NULL, shape = 19, size = 2, linetype = "solid", 
          bins = 30, 
          x_margin.plot = c("density", "histogram", "boxplot"), 
          x_margin.params = list(), 
          x_margin.ggtheme = theme_void(), 
          x_margin.space = FALSE, 
          x_margin.plot.size = 1, 
          
          plot_xmargin = TRUE,

          y_margin.plot = c("density", "histogram", "boxplot"), 
          y_margin.params = list(), 
          y_margin.ggtheme = theme_void(), 
          y_margin.space = FALSE, 
          y_margin.plot.size = 1, 
          
          plot_ymargin = TRUE,
          
          ymargin_as_xmargin = FALSE,
                    
          main.plot.size = 2, 
          title = NULL, xlab = NULL, 
          ylab = NULL, legend = "top", ggtheme = theme_pubr(),
          
          x_tit_font = c(12, "plain", "black"),
          x_text_font = c(12, "plain", "black"),
          y_tit_font = c(12, "plain", "black"),
          y_text_font = c(12, "plain", "black"),
          
          
          global_margins_cm = c(0, 0, 0.25, 0.25), # hard-coded default in ggscatterhist
          
          ...) 
{
  
  has_cowplot_v0.9 <- ggpubr:::has_cowplot_v0.9
  .check_margin_params <- ggpubr:::.check_margin_params
  .add_item <- ggpubr:::.add_item
  .insert_xaxis_grob <- ggpubr:::.insert_xaxis_grob
  .insert_yaxis_grob <- ggpubr:::.insert_yaxis_grob
  mutate <- dplyr:::mutate
  
  if (!has_cowplot_v0.9()) {
    warning("Install the latest developmental version of cowplot on github ", 
            "to fully use all the feature of ggscatterhist", 
            .call = FALSE)
  }
  x_margin.plot <- match.arg(x_margin.plot)
  y_margin.plot <- match.arg(y_margin.plot)
  
  x_margin.params <- .check_margin_params(x_margin.params, data, 
                                        color, fill, linetype)
  
  y_margin.params <- .check_margin_params(y_margin.params, data, 
                                          color, fill, linetype)
  
  
  if (x_margin.plot == "histogram") {
    x_margin.params <- x_margin.params %>% .add_item(bins = bins, 
                                                 position = "identity")
  }
  
  if (y_margin.plot == "histogram") {
    y_margin.params <- y_margin.params %>% .add_item(bins = bins, 
                                                     position = "identity")
  }
  
  if (!is.null(group)) {
    if (missing(color)) 
      color = group
    if (missing(shape)) 
      shape = group
  }
  . <- NULL
  sp <- ggscatter(data, x, y, color = color, fill = fill, palette = palette, 
                  shape = shape, size = size, xlab = xlab, ylab = ylab, 
                  ggtheme = ggtheme, title = title, legend = legend, ...)
  
  sp <- ggpar(sp,
              font.x = x_tit_font,
              font.y = y_tit_font,
              font.xtickslab = x_text_font,
              font.ytickslab = y_text_font
              )
  
  
  x_geomfunc <- switch(x_margin.plot, histogram = geom_histogram, 
                     density = geom_density, boxplot = geom_boxplot, geom_histogram)
  
  y_geomfunc <- switch(y_margin.plot, histogram = geom_histogram, 
                       density = geom_density, boxplot = geom_boxplot, geom_histogram)
  
  if (x_margin.plot %in% c("density", "histogram")) {
    x_xplot.x <- x
    x_xplot.y <- NULL
    x_yplot.x <- y
    x_yplot.y <- NULL
  }   else if (x_margin.plot %in% c("boxplot")) {
    if (is.null(group)) {
      data <- data %>% mutate(.xgroupx. = factor(1))
      group = ".xgroupx."
    }
    x_xplot.x <- group
    x_xplot.y <- x
    x_yplot.x <- group
    x_yplot.y <- y
  }
  
  if (y_margin.plot %in% c("density", "histogram")) {
    y_xplot.x <- x
    y_xplot.y <- NULL
    y_yplot.x <- y
    y_yplot.y <- NULL
  }  else if (y_margin.plot %in% c("boxplot")) {
    if (is.null(group)) {
      data <- data %>% mutate(.xgroupx. = factor(1))
      group = ".xgroupx."
    }
    y_xplot.x <- group
    y_xplot.y <- x
    y_yplot.x <- group
    y_yplot.y <- y
  }
  
  
  x_xplot <- ggplot() + x_margin.params %>% .add_item(geomfunc = x_geomfunc, 
                                                  data = data, x = x_xplot.x, y = x_xplot.y, alpha = 0.7) %>% 
    do.call(geom_exec, .)
  x_xplot <- set_palette(x_xplot, palette)
  x_yplot <- ggplot() + x_margin.params %>% .add_item(geomfunc = x_geomfunc, 
                                                  data = data, x = x_yplot.x, y = x_yplot.y, alpha = 0.7) %>% 
    do.call(geom_exec, .)
  x_yplot <- set_palette(x_yplot, palette)
  if (x_margin.plot %in% c("density", "histogram")) 
    x_yplot <- x_yplot + coord_flip()
  else if (x_margin.plot %in% c("boxplot")) 
    x_xplot <- x_xplot + coord_flip()
  
  y_xplot <- ggplot() + y_margin.params %>% .add_item(geomfunc = y_geomfunc, 
                                                      data = data, x = y_xplot.x, y = y_xplot.y, alpha = 0.7) %>% 
    do.call(geom_exec, .)
  y_xplot <- set_palette(y_xplot, palette)
  y_yplot <- ggplot() + y_margin.params %>% .add_item(geomfunc = y_geomfunc, 
                                                      data = data, x = y_yplot.x, y = y_yplot.y, alpha = 0.7) %>% 
    do.call(geom_exec, .)
  y_yplot <- set_palette(y_yplot, palette)
  if (y_margin.plot %in% c("density", "histogram")) 
    y_yplot <- y_yplot + coord_flip()
  else if (y_margin.plot %in% c("boxplot")) 
    y_xplot <- y_xplot + coord_flip()
  
  .legend <- get_legend(sp)
  # sp <- sp + theme(plot.margin = grid::unit(c(0, 0, 0.25, 0.25), 
  #                                           "cm"))
  
  sp <- sp + theme(plot.margin = grid::unit(global_margins_cm, 
                                            "cm"))
  
  
  
  
  x_xplot <- x_xplot + x_margin.ggtheme + clean_theme() + rremove("legend") + 
    theme(plot.margin = grid::unit(c(0, 0, 0, 0), "cm"))
  x_yplot <- x_yplot + x_margin.ggtheme + clean_theme() + rremove("legend") + 
    theme(plot.margin = grid::unit(c(0, 0, 0, 0), "cm"))
  
  y_xplot <- y_xplot + y_margin.ggtheme + clean_theme() + rremove("legend") + 
    theme(plot.margin = grid::unit(c(0, 0, 0, 0), "cm"))
  y_yplot <- y_yplot + y_margin.ggtheme + clean_theme() + rremove("legend") + 
    theme(plot.margin = grid::unit(c(0, 0, 0, 0), "cm"))
  
  
  ##TMP
  margin.space <- x_margin.space
  
  final_xplot <- x_xplot
  rm(x_xplot)
  rm(y_xplot)
  
  if( ymargin_as_xmargin) {
    final_yplot <- x_yplot
  } else {
    final_yplot <- y_yplot
    y_margin.plot.size <- x_margin.plot.size 
  }
  
  rm(y_yplot)
  rm(x_yplot)
  
  if (margin.space) {
    common.legend <- FALSE
    if (!is.null(.legend)) 
      common.legend = TRUE
    
    sp <- sp + theme(plot.title = element_blank(), plot.subtitle = element_blank())
    
    if(plot_xmargin & plot_ymargin) {
      # arrange with both x and y margin plots
      fig <- ggarrange(final_xplot, NULL, sp, final_yplot, ncol = 2, nrow = 2, 
                       align = "hv", 
                       widths = c(main.plot.size, y_margin.plot.size), 
                       heights = c(x_margin.plot.size, main.plot.size), 
                       common.legend = common.legend, 
                       legend = legend)
      
    } else  if(plot_xmargin & ! plot_ymargin) {
      
      # arrange with only x margin plots
      # fig <- ggarrange(final_xplot, NULL, sp, NULL, ncol = 2, nrow = 2, 
      #                  align = "hv", 
      #                  widths = c(main.plot.size, 0), 
      #                  heights = c(x_margin.plot.size, main.plot.size), 
      #                  common.legend = common.legend, 
      #                  legend = legend)
      fig <- ggarrange(final_xplot, sp, ncol = 1, nrow = 2, 
                       align = "hv", 
                       widths = c(main.plot.size), 
                       heights = c(x_margin.plot.size, main.plot.size), 
                       common.legend = common.legend, 
                       legend = legend)
      
    } else  if(! plot_xmargin & plot_ymargin) {
      
      # arrange with only y margin plots
      # fig <- ggarrange(NULL, NULL, sp, final_xplot, ncol = 2, nrow = 2,
      #                  align = "hv",
      #                  widths = c(main.plot.size, y_margin.plot.size), 
      #                  heights = c(0, main.plot.size), 
      #                  common.legend = common.legend,
      #                  legend = legend)
      
      fig <- ggarrange(sp, final_yplot, ncol = 2, nrow = 1, 
                       align = "hv", 
                       widths = c(main.plot.size, y_margin.plot.size), 
                       heights = c(main.plot.size), 
                       common.legend = common.legend, 
                       legend = legend)
      
    } else  if(! plot_xmargin & plot_ymargin) {
      fig <- ggarrange(sp, ncol = 1, nrow = 1, 
                       align = "hv", 
                       widths = c(main.plot.size),
                       heights = c(main.plot.size), 
                       common.legend = common.legend, 
                       legend = legend)
    }
    if (!is.null(title)) 
      fig <- annotate_figure(fig, top = text_grob(title, 
                                                  color = "black", size = 13, face = "bold"))
  } 	 else {
    
    if(plot_xmargin) {
      fig <- .insert_xaxis_grob(sp, final_xplot, grid::unit(x_margin.plot.size/5, 
                                                      "null"), position = "top")
    } else {
      fig <- sp
    }
    
    if(plot_ymargin) {
      fig <- .insert_yaxis_grob(fig, final_yplot, grid::unit(y_margin.plot.size/5, 
                                                       "null"), position = "right")
      
    }  
    fig <- cowplot::ggdraw(fig)
  }
  fig
}

##################################################################
# in this version you can specify which variables you want to use as x- or y-axes in the marginal plots
my_ggscatterhist_v2 <- function (data, x, y, group = NULL, color = "black", fill = NA, 
                                 palette = NULL, shape = 19, size = 2, linetype = "solid", 
                                 bins = 30, 
                                 x_margin.plot = c("density", "histogram", "boxplot"), 
                                 x_margin.params = list(), 
                                 x_margin.ggtheme = theme_void(), 
                                 x_margin.space = FALSE, 
                                 x_margin.plot.size = 1, 
                                 
                                 my_x_xplot.x = NULL,
                                 my_x_xplot.y = NULL,
                                 my_x_yplot.x = NULL,
                                 my_x_yplot.y = NULL,
                                 
                                 plot_xmargin = TRUE,
                                 
                                 y_margin.plot = c("density", "histogram", "boxplot"), 
                                 y_margin.params = list(), 
                                 y_margin.ggtheme = theme_void(), 
                                 y_margin.space = FALSE, 
                                 y_margin.plot.size = 1, 
                                 
                                 my_y_xplot.x = NULL,
                                 my_y_xplot.y = NULL,
                                 my_y_yplot.x = NULL,
                                 my_y_yplot.y = NULL,
                                 
                                 plot_ymargin = TRUE,
                                 
                                 ymargin_as_xmargin = FALSE,
                                 
                                 main.plot.size = 2, 
                                 title = NULL, xlab = NULL, 
                                 ylab = NULL, legend = "top", ggtheme = theme_pubr(), 
                                 
                                 global_margins_cm = c(0, 0, 0.25, 0.25), # hard-coded default in ggscatterhist
                                 
                                 ...) 
{
  
  has_cowplot_v0.9 <- ggpubr:::has_cowplot_v0.9
  .check_margin_params <- ggpubr:::.check_margin_params
  .add_item <- ggpubr:::.add_item
  .insert_xaxis_grob <- ggpubr:::.insert_xaxis_grob
  .insert_yaxis_grob <- ggpubr:::.insert_yaxis_grob
  mutate <- dplyr:::mutate
  
  if (!has_cowplot_v0.9()) {
    warning("Install the latest developmental version of cowplot on github ", 
            "to fully use all the feature of ggscatterhist", 
            .call = FALSE)
  }
  x_margin.plot <- match.arg(x_margin.plot)
  y_margin.plot <- match.arg(y_margin.plot)
  
  x_margin.params <- .check_margin_params(x_margin.params, data, 
                                          color, fill, linetype)
  
  y_margin.params <- .check_margin_params(y_margin.params, data, 
                                          color, fill, linetype)
  
  if (x_margin.plot == "histogram") {
    x_margin.params <- x_margin.params %>% .add_item(bins = bins, 
                                                     position = "identity")
  }
  if (y_margin.plot == "histogram") {
    y_margin.params <- y_margin.params %>% .add_item(bins = bins, 
                                                     position = "identity")
  }
  
  if (!is.null(group)) {
    if (missing(color)) 
      color = group
    if (missing(shape)) 
      shape = group
  }
  . <- NULL
  sp <- ggscatter(data, x, y, color = color, fill = fill, palette = palette, 
                  shape = shape, size = size, xlab = xlab, ylab = ylab, 
                  ggtheme = ggtheme, title = title, legend = legend, ...)
  
  x_geomfunc <- switch(x_margin.plot, histogram = geom_histogram, 
                       density = geom_density, boxplot = geom_boxplot, geom_histogram)
  
  y_geomfunc <- switch(y_margin.plot, histogram = geom_histogram, 
                       density = geom_density, boxplot = geom_boxplot, geom_histogram)
  
  if (x_margin.plot %in% c("density", "histogram")) {
    # x_xplot.x <- x
    x_xplot.x <- ifelse(is.null(my_x_xplot.x), x, my_x_xplot.x)
    x_xplot.y <- NULL
    # x_yplot.x <- y
    x_yplot.x <- ifelse(is.null(my_x_yplot.x), y, my_x_yplot.x) 
    x_yplot.y <- NULL
  }   else if (x_margin.plot %in% c("boxplot")) {
    if (is.null(group)) {
      data <- data %>% mutate(.xgroupx. = factor(1))
      group = ".xgroupx."
    }
    x_xplot.x <- group
    # x_xplot.y <- x
    x_xplot.y <- ifelse(is.null(my_x_xplot.y), x, my_x_xplot.y)
    x_yplot.x <- group
    # x_yplot.y <- y
    x_yplot.y <- ifelse(is.null(my_x_yplot.y), y, my_x_yplot.y)
  }
  
  if (y_margin.plot %in% c("density", "histogram")) {
    # y_xplot.x <- x
    y_xplot.x <- ifelse(is.null(my_y_xplot.x), x, my_y_xplot.x)
    y_xplot.y <- NULL
    # y_yplot.x <- y
    y_yplot.x <- ifelse(is.null(my_y_yplot.x), y, my_y_yplot.x)
    y_yplot.y <- NULL
  }  else if (y_margin.plot %in% c("boxplot")) {
    if (is.null(group)) {
      data <- data %>% mutate(.xgroupx. = factor(1))
      group = ".xgroupx."
    }
    y_xplot.x <- group
    # y_xplot.y <- x
    y_xplot.y <- ifelse(is.null(my_y_xplot.y), x, my_y_xplot.y)
    y_yplot.x <- group
    # y_yplot.y <- y
    y_yplot.y <- ifelse(is.null(my_y_yplot.y), y, my_y_yplot.y)
  }
  
  x_xplot <- ggplot() + x_margin.params %>% .add_item(geomfunc = x_geomfunc, 
                                                      data = data, x = x_xplot.x, y = x_xplot.y, alpha = 0.7) %>% 
    do.call(geom_exec, .)
  x_xplot <- set_palette(x_xplot, palette)
  x_yplot <- ggplot() + x_margin.params %>% .add_item(geomfunc = x_geomfunc, 
                                                      data = data, x = x_yplot.x, y = x_yplot.y, alpha = 0.7) %>% 
    do.call(geom_exec, .)
  x_yplot <- set_palette(x_yplot, palette)
  if (x_margin.plot %in% c("density", "histogram")) 
    x_yplot <- x_yplot + coord_flip()
  else if (x_margin.plot %in% c("boxplot")) 
    x_xplot <- x_xplot + coord_flip()
  
  
  y_xplot <- ggplot() + y_margin.params %>% .add_item(geomfunc = y_geomfunc, 
                                                      data = data, x = y_xplot.x, y = y_xplot.y, alpha = 0.7) %>% 
    do.call(geom_exec, .)
  y_xplot <- set_palette(y_xplot, palette)
  y_yplot <- ggplot() + y_margin.params %>% .add_item(geomfunc = y_geomfunc, 
                                                      data = data, x = y_yplot.x, y = y_yplot.y, alpha = 0.7) %>% 
    do.call(geom_exec, .)
  y_yplot <- set_palette(y_yplot, palette)
  if (y_margin.plot %in% c("density", "histogram")) 
    y_yplot <- y_yplot + coord_flip()
  else if (y_margin.plot %in% c("boxplot")) 
    y_xplot <- y_xplot + coord_flip()
  
  .legend <- get_legend(sp)
  # sp <- sp + theme(plot.margin = grid::unit(c(0, 0, 0.25, 0.25), 
  #                                           "cm"))
  sp <- sp + theme(plot.margin = grid::unit(global_margins_cm, 
                                            "cm"))
  
  x_xplot <- x_xplot + x_margin.ggtheme + clean_theme() + rremove("legend") + 
    theme(plot.margin = grid::unit(c(0, 0, 0, 0), "cm"))
  x_yplot <- x_yplot + x_margin.ggtheme + clean_theme() + rremove("legend") + 
    theme(plot.margin = grid::unit(c(0, 0, 0, 0), "cm"))
  
  y_xplot <- y_xplot + y_margin.ggtheme + clean_theme() + rremove("legend") + 
    theme(plot.margin = grid::unit(c(0, 0, 0, 0), "cm"))
  y_yplot <- y_yplot + y_margin.ggtheme + clean_theme() + rremove("legend") + 
    theme(plot.margin = grid::unit(c(0, 0, 0, 0), "cm"))
  
  ##TMP !!! margin.space use x_margin.space !!!
  margin.space <- x_margin.space
  
  final_xplot <- x_xplot
  rm(x_xplot)
  rm(y_xplot)
  
  
  if( ymargin_as_xmargin) {
    final_yplot <- x_yplot
  } else {
    final_yplot <- y_yplot
    y_margin.plot.size <- x_margin.plot.size 
  }
  
  rm(y_yplot)
  rm(x_yplot)
  
  if (margin.space) {
    common.legend <- FALSE
    if (!is.null(.legend)) 
      common.legend = TRUE
    
    sp <- sp + theme(plot.title = element_blank(), plot.subtitle = element_blank())
    
    if(plot_xmargin & plot_ymargin) {
      # arrange with both x and y margin plots
      fig <- ggarrange(final_xplot, NULL, sp, final_yplot, ncol = 2, nrow = 2, 
                       align = "hv", 
                       widths = c(main.plot.size, y_margin.plot.size), 
                       heights = c(x_margin.plot.size, main.plot.size), 
                       common.legend = common.legend, 
                       legend = legend)
      
    } else  if(plot_xmargin & ! plot_ymargin) {
      
      # arrange with only x margin plots
      # fig <- ggarrange(final_xplot, NULL, sp, NULL, ncol = 2, nrow = 2, 
      #                  align = "hv", 
      #                  widths = c(main.plot.size, 0), 
      #                  heights = c(x_margin.plot.size, main.plot.size), 
      #                  common.legend = common.legend, 
      #                  legend = legend)
      fig <- ggarrange(final_xplot, sp, ncol = 1, nrow = 2, 
                       align = "hv", 
                       widths = c(main.plot.size), 
                       heights = c(x_margin.plot.size, main.plot.size), 
                       common.legend = common.legend, 
                       legend = legend)
      
    } else  if(! plot_xmargin & plot_ymargin) {
      
      # arrange with only y margin plots
      # fig <- ggarrange(NULL, NULL, sp, final_xplot, ncol = 2, nrow = 2,
      #                  align = "hv",
      #                  widths = c(main.plot.size, y_margin.plot.size), 
      #                  heights = c(0, main.plot.size), 
      #                  common.legend = common.legend,
      #                  legend = legend)
      fig <- ggarrange(sp, final_yplot, ncol = 2, nrow = 1, 
                       align = "hv", 
                       widths = c(main.plot.size, y_margin.plot.size), 
                       heights = c(main.plot.size), 
                       common.legend = common.legend, 
                       legend = legend)
      
    } else  if(! plot_xmargin & plot_ymargin) {
      fig <- ggarrange(sp, ncol = 1, nrow = 1, 
                       align = "hv", 
                       widths = c(main.plot.size),
                       heights = c(main.plot.size), 
                       common.legend = common.legend, 
                       legend = legend)
    }
    
    if (!is.null(title)) 
      fig <- annotate_figure(fig, top = text_grob(title, 
                                                  color = "black", size = 13, face = "bold"))
  }  else {
    if(plot_xmargin) {
      fig <- .insert_xaxis_grob(sp, final_xplot, grid::unit(x_margin.plot.size/5, 
                                                            "null"), position = "top")
    } else {
      fig <- sp
    }
    if(plot_ymargin) {
      fig <- .insert_yaxis_grob(fig, final_yplot, grid::unit(y_margin.plot.size/5, 
                                                             "null"), position = "right")
    }  
    fig <- cowplot::ggdraw(fig)
  }

  fig
}


# my_ggscatterhist(
#   family_dist_coexpr_DT, 
#   x = "dist_kb", 
#   xlab=my_xlab,
#   ylab=my_ylab,
#   y = "coexpr",
#   point = FALSE,
#   # rug=TRUE,
#   add = "loess",
#   color = "sameTAD_lab", size = 3, alpha = 0.6,
#   palette = c("darkslateblue", "darkorange1"),
#   x_margin.params = list(fill = "sameTAD_lab", color = "black", size = 0.2),
#   x_margin.ggtheme = theme_minimal(),
#   x_margin.plot = "density",
#   y_margin.params = list(fill = "sameTAD_lab", color = "black", size = 0.2),
#   y_margin.ggtheme = theme_minimal(),
#   y_margin.plot = "boxplot",
#   plot_xmargin = TRUE,
#   plot_ymargin=TRUE,
#   ymargin_as_xmargin = FALSE
# )
# my_ggscatterhist_v2(
#   family_dist_coexpr_DT, 
#   x = "dist_kb", 
#   xlab=my_xlab,
#   ylab=my_ylab,
#   y = "coexpr",
#   point = FALSE,
#   # rug=TRUE,
#   add = "loess",
#   color = "sameTAD_lab", size = 3, alpha = 0.6,
#   palette = c("darkslateblue", "darkorange1"),
#   x_margin.params = list(fill = "sameTAD_lab", color = "black", size = 0.2),
#   x_margin.ggtheme = theme_minimal(),
#   x_margin.plot = "density",
#   y_margin.params = list(fill = "sameTAD_lab", color = "black", size = 0.2),
#   y_margin.ggtheme = theme_minimal(),
#   y_margin.plot = "boxplot",
#   plot_xmargin = TRUE,
#   plot_ymargin=TRUE,
#   ymargin_as_xmargin = FALSE
# )

# my_ggscatterhist_v2(
#   family_dist_coexpr_DT, 
#   x = "dist_kb", 
#   xlab=my_xlab,
#   ylab=my_ylab,
#   y = "coexpr",
#   point = FALSE,
#   # rug=TRUE,
#   add = "loess",
#   color = "sameTAD_lab", size = 3, alpha = 0.6,
#   palette = c("darkslateblue", "darkorange1"),
#   x_margin.params = list(fill = "sameTAD_lab", color = "black", size = 0.2),
#   x_margin.ggtheme = theme_minimal(),
#   x_margin.plot = "density",
#   y_margin.params = list(fill = "sameTAD_lab", color = "black", size = 0.2),
#   y_margin.ggtheme = theme_minimal(),
#   y_margin.plot = "boxplot",
#   plot_xmargin = TRUE,
#   plot_ymargin=TRUE,
#   ymargin_as_xmargin = FALSE,
#   my_y_yplot.y = "dist_kb", # reverse the intuitive reference
#   my_x_xplot.x = "coexpr"  # reverse the intuitive reference
# )

# group = NULL
# color = "black"
# fill = NA
# palette = NULL
# shape = 19
# size = 2
# linetype = "solid"
# bins = 30
# x_margin.plot = c("density", "histogram", "boxplot")
# x_margin.params = list()
# x_margin.ggtheme = theme_void()
# x_margin.space = FALSE
# x_margin.plot.size = 1 
# plot_xmargin = TRUE
# y_margin.plot = c("density", "histogram", "boxplot")
# y_margin.params = list()
# y_margin.ggtheme = theme_void()
# y_margin.space = FALSE
# y_margin.plot.size = 1 
# plot_ymargin = TRUE
# ymargin_as_xmargin = FALSE
# main.plot.size = 2
# title = NULL
# xlab = NULL
# ylab = NULL
# legend = "top"
# ggtheme = theme_pubr()
# 
# 
#   data=family_dist_coexpr_DT
#   x = "dist_kb"
#   xlab=my_xlab
#   ylab=my_ylab
#   title = paste0("Gene pair expr. corr. vs. dist. - ", curr_dataset)
#   subtitle=paste0("(", i_fam," - nTop=", nTopTADs, ")")
#   legend.title=""
#   y = "coexpr"
#   point = FALSE
#   # rug=TRUE,
#   add = fitMeth
#   color = "sameTAD_lab"
#   size = 3
#   alpha = 0.6
#   palette = c(col1, col2)
#   x_margin.params = list(fill = "sameTAD_lab", color = "black", size = 0.2)
#   x_margin.ggtheme = theme_minimal()
#   x_margin.plot = "density"
#   y_margin.params = list(fill = "sameTAD_lab", color = "black", size = 0.2)
#   y_margin.ggtheme = theme_minimal()
#   y_margin.plot = "boxplot"
#   plot_xmargin = TRUE
#   plot_ymargin=TRUE
#   ymargin_as_xmargin = FALSE

