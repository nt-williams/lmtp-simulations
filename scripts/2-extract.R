# Nick Williams
# Research Biostatistician
# Department of Population Health Sciences
# Weill Cornell Medicine

# packages ----------------------------------------------------------------

library(lmtp)
library(data.table)
library(purrr)

# global ------------------------------------------------------------------

conf <- config::get()

set.seed(48838)

seed <- sample(1e8, conf$reps)
task <- setDT(expand.grid(seed = seed, type = 1:4, n = conf$nobs))

res <- here::here("results")

read_files <- function(pattern) {
  files <- list.files(res, pattern = paste0("^", pattern))
  task <- gsub("(.+)_([0-9]+).rds", "\\2", files)
  out <- lapply(files, function(x) readRDS(file.path(res, x)))
  names(out) <- task
  return(out)
}

find_errors <- function(x) {
  er <- which(as.vector(unlist(lapply(x, inherits, what = "try-error"))))
  if (length(er) == 0) return(x)
  x[-er]
}

tidy_and_bind <- function(x) {
  map_dfr(x, tidy, .id = "task")
}

# data import -------------------------------------------------------------

tml <- read_files("tml")
sdr <- read_files("sdr")
sub <- read_files("sub")
ipw <- read_files("ipw")

fits <- list(tml = tml, 
             sdr = sdr, 
             sub = sub, 
             ipw = ipw)

# extract -----------------------------------------------------------------

fits <- map(fits, find_errors) %>% 
  map_dfr(tidy_and_bind) %>% 
  setDT()

task[, task := .I]
fits[, task := as.numeric(task)]
setkey(task, task)
setkey(fits, task)

fits <- task[fits]

# export ------------------------------------------------------------------

saveRDS(fits, file.path(res, "fits-exported.rds"))
