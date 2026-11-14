This repository provides codes for implementing the procedure in the paper "Detecting Multiple Change Points in Linear Models With Heteroscedasticity", by Lajos Horvath, Gregory Rice, Yuqian Zhao (2026, Journal of Time Series Analysis)

# DGP + V/Q Tests 
# `dgp.R`: Main simulation setting with Normal, AR, GARCH, or smoothly heteroscedastic errors. 
# `vq_test.R`: V/Q tests and 5% critical values.
# `demo.R`: Simple examples with and without a coefficient change.
## Main Setting  The initial coefficients are `(1,1)`. Both increase by `delta` after `floor(N*kstar)`; `delta=0` means no coefficient change. 
## Defaults: `N=250`, `kstar=0.5`.  The DGP follows the main setting in the original `dgp.R`. 
## Default V/Q results preserve the original implementation and have passed numerical equivalence checks.
