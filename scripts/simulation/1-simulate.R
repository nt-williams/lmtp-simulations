
# Nick Williams
# Research Biostatistician
# Department of Population Health Sciences
# Weill Cornell Medicine

.libPaths("/home/niw4001/R_local")

library(lmtp)
library(future)

conf <- config::get()

source(here::here("R", "generation.R"))
source(here::here("R", "SL.hal.R"))
source(here::here("R", "simulate.R"))

set.seed(48838)

seed <- sample(1e8, conf$reps)
task <- expand.grid(seed = seed, type = 1:4, n = conf$nobs)

plan(multisession)

partition(task, conf$task, conf$machines, conf$save, "tml")
partition(task, conf$task, conf$machines, conf$save, "sdr")
partition(task, conf$task, conf$machines, conf$save, "sub")
partition(task, conf$task, conf$machines, conf$save, "ipw")

quit("no")
