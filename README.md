# DGP-theory

The goal of the `dgp` package is computing the posterior density of local extrema shown in the paper ``Semiparametric Bayesian inference for local extrema of functions in the presence of noise".


## Installation and Usage

You can install `dgp` from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("chenghanyustats/dgp")
```

Some important functions are

- `Rsolnp::solnp()` is used to optimize the hyperparameters of the kernel function used in the GPR.

- `log_post_t_theory()` computes the log posterior density of local extrema.

- `get_pred_ci_gp()` does the curve fitting of DGP regression and computes the associated credible interval.

- `plot_pred_gp_f_y()` plots the fitted results from `get_pred_ci_gp()`.

- `get_hpd_interval_from_den()` gets the HPD interval from the posterior density of local extrema.

- `get_map()` gets the maximum a posteriori and estimates the number of stationary points.


## Package Dependency

The required R packages include **matrixcalc**, **Rsolnp**, **emulator**, and **mvnfast**. The package **doParallel** is used for simple multi-core computing. The **ftnonpar** and **KernSmooth** packages are required for the smoothed taut string methood (STS) and nonparametric kernel smoothing (NKS) method in the paper, respectively. The R version >= 3.5 is required.



### Simulation Study

The package provide a vignette `simulation.Rmd`(`simulation.pdf`) for implementing the proposed method and conducting a simulation study. The simulated data are saved in `./data` that are generated from `./data-raw/gen_sim_data.R`.



### ERP Data Analysis

The package provide a vignette `erpdata.Rmd`(`erpdata.pdf`) for implementing the proposed method and conducting a simulation study. The raw ERP data is saved in `./data/Raw_ERP.csv`.

