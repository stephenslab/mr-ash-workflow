library(dscrutils)
dsc <- dscquery("mr_ash",
                targets = c("simulate","simulate.s","simulate.sigma",
                            "fit","fit.lambda_est_method"))
dsc <- transform(dsc,
                 simulate = factor(simulate),
                 fit      = factor(fit),
                 fit.lambda_est_method = factor(fit.lambda_est_method))
summary(dsc[c("simulate","fit","fit.lambda_est_method")])
