
library(tidyverse)
library(patchwork)
library(here)
library(colorspace)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# load data, run bootstrap, load functions ----
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

rm(list=ls()) # clear workspace
source("R/01_import_data.R")
source("R/02_bootstrap.R")
source("R/03_functions.R")

#@@@@@@@@@@@@@@@@@@@@@@
# Set parameters ----
#@@@@@@@@@@@@@@@@@@@@@

# set upper and lower bounds of otolith size to integrate over
lower <- 0
upper <- max(useData$OL.t1, na.rm = TRUE)+0.2
n_breaks <- 5000

# set max age
max_age <- max(useData$age)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Simple mean projection -----
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


proj_mean <-
  projectO(
    iparsO = c(l = lower, u = upper, n = n_breaks), #set upper and lower bounds of otolith size
    mpars = apply(mpars_bs, 2, mean),
    ref_temp = mean(unique(useData$covariate)),
    max_a = max_age
  )

## plot the expected age specific density functions
par(mfrow = c(1, 1))
with(proj_mean, {
  plot(
    mesh,
    n[, 1],
    type = "l",
    ylim = c(0, 15),
    xlab = "Otolith Size",
    ylab = "Density"
  )
  for (A in 2:max_a) lines(mesh, n[, A])
})

## plot the size vs age relationship...
par(mfrow = c(1, 1))
with(proj_mean, {
  mean_sz <- apply(n, 2, function(pdens) sum(pdens * mesh * dlta))
  plot(
    1:max_a - 1,
    mean_sz,
    xlab = "Age",
    ylab = "Otolith Size",
    ylim = c(0, 0.7)
  )
  ## ...and add the 90% quantiles
  quantile90 <- apply(n, 2, get_quantile, dlta, mesh)
  for (A in 1:max_a - 1)
    lines(rep(A, length(quantile90[[A + 1]])), quantile90[[A + 1]])
})



## repeat the mean calculations for different temperatures (+/- 1 degree)
proj_mean_m1 <-
  projectO(
    iparsO = c(l = 0, u = 0.7, n = 500),
    mpars = apply(mpars_bs, 2, mean),
    ref_temp = mean(unique(useData$covariate)) - 1,
    max_a = 19
  )
proj_mean_p1 <-
  projectO(
    iparsO = c(l = 0, u = 0.7, n = 500),
    mpars = apply(mpars_bs, 2, mean),
    ref_temp = mean(unique(useData$covariate)) + 1,
    max_a = 19
  )

## plot the size vs age relationship...
par(mfrow = c(1, 1))
with(proj_mean, {
  mean_sz <- apply(n, 2, function(pdens) sum(pdens * mesh * dlta))
  plot(
    1:max_a - 1,
    mean_sz,
    ylim = c(0, 0.7),
    xlab = "Age",
    ylab = "Otolith Size",
    type = "b"
  )
})
with(proj_mean_m1, {
  mean_sz <- apply(n, 2, function(pdens) sum(pdens * mesh * dlta))
  lines(1:max_a - 1, mean_sz, type = "b", col = "blue")
})
with(proj_mean_p1, {
  mean_sz <- apply(n, 2, function(pdens) sum(pdens * mesh * dlta))
  lines(1:max_a - 1, mean_sz, type = "b", col = "red")
})