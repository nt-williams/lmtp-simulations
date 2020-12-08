## LMTP Simulations

> Simulation material for ["Non-parametric causal effects based on longitudinal modified treatment policies."](https://arxiv.org/abs/2006.01366)

**Authors:** Ivan Díaz, Nicholas Williams, Katherine Hoffman, and Edward Schenck

------------------------------------------------------------------------

### Description

This repository contains the code used for the simulation study in the paper "Non-parametric casual effects based on longitudinal modified treatment policies." Simulations are designed to be run on a high performance computing cluster using Slurm. The main file structure is as follows:

-   R: Contains R function definition scripts.

-   scripts: Contains R scripts that run the simulations and summarize simulation results.

-   results: Contains simulation results.

**This simulation uses a special version of the [*lmtp*](https://github.com/nt-williams/lmtp)** **package that allows for different [*SuperLearner*](https://github.com/ecpolley/SuperLearner) calls at every time point.** This package version is available in this repository as "lmtp\_1.0.0.9001.tar.gz". It can be installed using:

```{r, eval = FALSE}
install.packages("lmtp_1.0.0.9001.tar.gz", repos = NULL, type = "source")
```
