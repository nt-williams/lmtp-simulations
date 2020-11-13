
# Nick Williams
# Research Biostatistician 
# Department of Healthcare Policy and Research 
# Weill Cornell Medicine

.libPaths("/home/niw4001/R_local")

conf <- config::get(config = "truth")

sources <- here::here("R", list.files(here::here("R")))
sapply(sources, source, .GlobalEnv)

set.seed(5253)
seed <- sample(1e5, conf$reps)[conf$instance]

partition_truth(conf$nobs, conf$instance, seed, conf$save)

quit("no")
