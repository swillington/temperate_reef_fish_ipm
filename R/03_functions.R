#@@@@@@@@@@@@@@@@@@@@@@@
## Projection functions
#@@@@@@@@@@@@@@@@@@@@@@@

# size kernel function
szkern <- function(TL, z, pars) {
  TLmean <- pars["smod.(Intercept)"] + pars["smod.rad4"] * z
  TLsd <- pars["smod.sd.p0"] * sqrt(exp(2 * pars["smod.sd.p"] * z))
  return(dnorm(TL, mean = TLmean, sd = TLsd)) # returns probability density function for body size
}

# growth kernel function
gkern <- function(z1, z, T, pars) {
  z1mean <- pars["gmod.(Intercept)"] +
    pars["gmod.OL.t"] * z +
    pars["gmod.covariate"] * T
  z1sd <- pars["gmod.sd"]
  return(dnorm(z1, mean = z1mean, sd = z1sd))
}

# initial annuli/increment - set the baseline. Predicted based on temperature only, as cannot predict using the previous increment
idist <- function(z1, T, pars) {
  z1mean <- pars["imod.(Intercept)"] + pars["imod.covariate"] * T
  z1sd <- pars["imod.sd"]
  return(dnorm(z1, mean = z1mean, sd = z1sd))
}

# width to integrate over
mkdelta <- function(ipars) {
  (ipars["u"] - ipars["l"]) / ipars["n"]
}

# create mesh to integrate over
mkmesh <- function(ipars) {
  ipars["l"] + mkdelta(ipars) * (1:ipars["n"] - 1 / 2)
}

# getting the 95th quantile of each probability density function
get_quantile <- function(pdens, dlta, mesh, quantile = 0.95) {
  cdf <- cumsum(pdens) / sum(pdens)
  range(mesh[cdf > (1 - quantile) / 2 & cdf < 1 / 2 + quantile / 2])
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# otolith projection function
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@

projectO <- function(iparsO, mpars, max_a, ref_temp) {
  mesh <- mkmesh(iparsO)
  dlta <- mkdelta(iparsO)
  gK <- outer(mesh, mesh, gkern, T = ref_temp, pars = mpars) * dlta
  n <- matrix(NA, nrow = length(mesh), ncol = max_a)
  n[, 1] <- idist(mesh, T = ref_temp, pars = mpars)
  for (A in 2:max_a) {
    n[, A] <- (gK %*% n[, A - 1])[,, drop = TRUE]
  }
  colnames(n) <- paste0("age_", 1:max_a - 1)
  return(list(
    n = n, # <- projections
    ref_temp = ref_temp, # <- environment
    max_a = max_a,
    mesh = mesh,
    dlta = dlta
  )) # <- integr. params
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# body size projection function 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

projectB <- function(iparsO, iparsB, mpars, max_a, ref_temp) {
  meshB <- mkmesh(iparsB)
  dltaB <- mkdelta(iparsB)
  meshO <- mkmesh(iparsO)
  dltaO <- mkdelta(iparsO)
  
  # making the growth (# 1) and size (# 2) kernels
  gK <- outer(meshO, meshO, gkern, T = ref_temp, pars = mpars) * dltaO # 1
  bK <- outer(meshB, meshO, szkern, pars = mpars) * dltaO # 2
  
  # creating empty matrix to put output into
  nO <- matrix(NA, nrow = length(meshO), ncol = max_a)
  nB <- matrix(NA, nrow = length(meshB), ncol = max_a)
  
  # otolith size (nO) and fish body size (nB) distributions of age0 individuals:
  nO[, 1] <- idist(meshO, T = ref_temp, pars = mpars)
  nB[, 1] <- (bK %*% nO[, 1])[,, drop = TRUE]
  
  # extend this distribution to all ages(using growth kernel and body kernel)
  for (A in 2:max_a) {
    nO[, A] <- (gK %*% nO[, A - 1])[,, drop = TRUE]
    nB[, A] <- (bK %*% nO[, A])[,, drop = TRUE]
  }
  colnames(nB) <- paste0("age_", 1:max_a - 1)
  colnames(nO) <- paste0("age_", 1:max_a - 1)
  return(list(
    nB = nB, # <- projections
    ref_temp = ref_temp, # <- environment
    max_a = max_a,
    meshB = meshB,
    dltaB = dltaB
  )) # <- integr. params
}