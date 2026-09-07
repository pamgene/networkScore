# networkScore

Computes golden score (paired kinase+sensitivity) and kinase-only score for
a condition, via permutation testing against networks built by
[`networkGen`](../networkGen). See the suite-wide
[Context Map](../networkGen/CONTEXT-MAP.md).

## Install

```r
# install.packages("remotes")
remotes::install_local("path/to/networkGen") # dependency; needs the real PCSF package too, see its README
remotes::install_local("path/to/networkScore")
```

## Getting started

See `vignette("networkScore")` for a full walkthrough on toy data.

```r
library(networkScore)

result <- make_golden_score(
  uka = uka, spec_cutoff = 1.0, perc_cutoffs = c(0.5, 0.7),
  ppi_network = ppi_network, b = 1.5, respath = "results/my_run"
)

result$results
```

Set a `future` plan before calling for actual parallelism across the batch
of networks a run builds (one observed + N permuted per condition):

```r
future::plan(future::multisession, workers = 4)
```
