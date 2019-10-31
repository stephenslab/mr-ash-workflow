library(dscrutils)
library(ggplot2)
library(cowplot)

# Load the DSC results.
dsc <- dscquery("mr_ash",
                targets = c("simulate","simulate.s","simulate.sigma",
                            "fit","fit.lambda_est_method","score.err"))
rows <- which(is.na(dsc$fit.lambda_est_method))
dsc$fit.lambda_est_method[rows] <- ""
dsc <- transform(dsc,fit = paste(fit,fit.lambda_est_method))
dsc <- transform(dsc,
                 simulate              = factor(simulate),
                 simualte.s            = factor(simulate.s),
                 fit                   = factor(fit),
                 fit.lambda_est_method = factor(fit.lambda_est_method))

# Summarize the prediction error separately for each setting of "s", the
# number of simulated effects.
p <- ggplot(dsc,aes(x = fit,y = score.err,fill = fit)) +
     geom_boxplot(color = "black",outlier.size = 1,width = 0.85) +
     facet_wrap(~simulate.s,nrow = 2,scales = "free") +
     labs(x = "method",y = "rsse") +
     theme_cowplot(font_size = 10) +
     theme(axis.text.x = element_text(angle = 45,hjust = 1))
