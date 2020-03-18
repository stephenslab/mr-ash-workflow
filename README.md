# MR.ASH (github: mr-ash-workflow)

Numerical experiments for evaluating the "multiple regression with adaptive
shrinkage" (MR.ASH) method and other linear regression methods that
are well suited for large-scale data sets.

This workflow is for reproducing the results used in the manuscript

[MR.ASH: A Novel Variational Empirical Bayes Approach to Multiple Linear Regression][mr-ash-manu].

---------------------------------------------------------------

## Reproducing the result using R

The proposed VEB (Variational Empirical Bayes) approach can be implemented via our companion package `mr.ash.alpha`,
which can be found in [MR.ASH github][mr-ash-alpha].

### Quick Start for mr.ash.alpha

To install the latest version of the R package `mr.ash.alpha` from [MR.ASH github][mr-ash-alpha], clone or
download the git repository, then use the `install_local` function from
`devtools`, replacing "path/to/repo" with the appropriate directory path:

```r
devtools::install_local("path/to/repo")
```

Also, one can install the package via

```r
devtools::install_github("stephenslab/mr.ash.alpha")
```

### Quick Start for Other comparision methods

The list of competitors is as follows.

```r
# Penalized linear regression based on Ridge, Lasso, Elastic Net penalties, fitted by CV (Ridge, Lass, E-NET)
install.packages("glmnet")
# Penalized linear regression based on SCAD, MCP penalties, fitted by CV (SCAD, MCP)
install.packages("ncvreg")
# Penalized linear regression based on L0 penalty (or L0+L1, L0+L2), fitted by CV (L0Learn)
install.packages("L0Learn")
# Bayesian linear regression based on additive single effect model, fitted by VB (SuSiE)
install_github("stephenslab/susieR")
# Bayesian linear regression based on hieararchical spike-and-slab prior, fitted by MCMC (BayesB, BLasso)
install.packages("BGLR")
# Bayesian linear regression based on hieararchical spike-and-slab prior, fitted by VB and discrete BMA (varbvs)
install.packages("varbvs")
```

For detailed discussion on the comparison methods, please see the package documentations and reference therein.

### Implementation of the methods

For each method, its implementation has a wrapper in `code/method_wrapper.R`.

For each simulation setting, its implementation has a wrapper in `code/sim_wrapper.R`.

### Results

We provide `*.Rmd` files for reproducing all the simulations in `analysis` folder.

We provide `*.RDS` files for storing all the results from `analysis/*.Rmd` in `results` folder.

We provide `analysis/plots_for_the_paper.RDS` for reproducing all the figures in our manuscript [MR.ASH][mr-ash-manu].

### Figures

| Figure     | Description                                                              | Source code for reproducing the result | Source code for plotting |   |
|------------|--------------------------------------------------------------------------|----------------------------------------|--------------------------|---|
| `Figure 1` | Flexibility of MR.ASH prior, in terms of shrinkage operator and penalty  | `analysis/Figure1_Sourcecode`          | `analysis/Figure1_Plot`  |   |
| `Figure 2` | Adaptation to Sparsity (Penalized Linear Regression)                     | `analysis/Figure2_Sourcecode`          | `analysis/Figure2_Plot`  |   |
| `Figure 3` | Sample size, signal strength, signal shape (Penalized Linear Regression) | `analysis/Figure3_Sourcecode`          | `analysis/Figure3_Plot`  |   |
| `Figure 4` | Adaptation to Sparsity (Bayesian Linear Regression)                      | `analysis/Figure2_Sourcecode`          | `analysis/Figure4_Plot`  |   |
| `Figure 5` | Sample size, signal strength, signal shape (Bayesian Linear Regression)  | `analysis/Figure3_Sourcecode`          | `analysis/Figure5_Plot`  |   |
| `Figure 6` | Computation time (scalability)                                           | `analysis/Figure6_Sourcecode`          | `analysis/Figure6_Plot`  |   |
| `Figure 7` | Relative Performance                                                     | -          | `analysis/Figure7_Plot`  |   |


---------------------------------------------------------------

## Reproducing the results using DSC

Not ready yet

---------------------------------------------------------------

## License

Copyright (c) 2019-2020, Youngseok Kim.

All source code and software in this repository are made available
under the terms of the [MIT license][mit-license]. See
file [LICENSE](LICENSE) for the full text of the license.

---------------------------------------------------------------

## Credits

The mr.ash R package was developed by [Youngseok Kim][youngseok].

[mit-license]: https://opensource.org/licenses/mit-license.html
[devtools]: https://github.com/r-lib/devtools
[uchicago]: https://www.uchicago.edu
[youngseok]: https://github.com/youngseok-kim
[mr-ash-manu]: https://https://stephenslab.uchicago.edu/
[mr-ash-alpha]: https://github.com/stephenslab/mr.ash.alpha
