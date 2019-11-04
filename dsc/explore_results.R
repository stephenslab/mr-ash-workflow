# Summarize the prediction error separately for each setting of "s", the
# number of simulated effects.
p <- ggplot(dsc,aes(x = fit,y = score.err,fill = fit)) +
     geom_boxplot(color = "black",outlier.size = 1,width = 0.85) +
     facet_wrap(~simulate.s,nrow = 2,scales = "free") +
     labs(x = "method",y = "rsse") +
     theme_cowplot(font_size = 10) +
     theme(axis.text.x = element_text(angle = 45,hjust = 1))
