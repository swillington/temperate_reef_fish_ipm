
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# set bootstrap to estimate error around parameters coeffs @
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


library(tidyverse)
library(nlme)
library(lme4)

# initial data
in_dat <-
  gr_dat %>%
  filter(annuli_num == 2) %>%
  select(FishID, OL.t, covariate) %>%
  drop_na()

# growth data
gr_dat <-
  gr_dat %>%
  select(FishID, OL.t1, OL.t, covariate) %>%
  filter(FishID %in% al_dat$FishID) %>% 
  drop_na()

# allometry data
al_dat <-
  al_dat %>%
  select(FishID, fishlength, rad4, age, yearcap) %>%
  drop_na()

#@@@@@@@@@@@@@@@@@@@@@
# set models here ####
#@@@@@@@@@@@@@@@@@@@@@

# initial model for first year of growth
im <- lm(OL.t ~ covariate, data = in_dat)

# growth model
gm <- lm(OL.t1 ~ OL.t + covariate, data = gr_dat)
gm2 <- lm(OL.t1 ~ OL.t*covariate, data = gr_dat)

# allometry model
am <- gls(fishlength ~ rad4, weights = varExp(form = ~rad4), data = al_dat)


summary(im)
summary(gm2)

#@@@@@@@@@@@@@@@@@@@@@@@@@
# set number of boots ####
#@@@@@@@@@@@@@@@@@@@@@@@@@

n_boots <- 1000


# empty list to store model estimates
in_ests <- gr_ests <- al_ests <- list()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# bootstrap loop to estimate error around model coefs
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

for (i in 1:n_boots) {
  # resample the growth data + fit and store the init size model coefs
  in_dat_now <- sample_n(in_dat, size = nrow(in_dat), replace = TRUE)
  # simple linear model of otolith increment width ~ temperature
  im <- im
  #store model coefs
  in_ests[[i]] <- as.list(c(imod = coef(im), imod.sd = summary(im)$sigma))
  
  # resample the growth data + fit and store the growth model coefs
  gr_dat_now <- sample_n(gr_dat, size = nrow(gr_dat), replace = TRUE)
  # linear model otolith increment at t+1 ~ otolith increment at t + temperature
  gm <- gm
  #store model coefs
  gr_ests[[i]] <- as.list(c(gmod = coef(gm), gmod.sd = summary(gm)$sigma))
  
  # fit and store the allometry model coefficients
  al_dat_now <- sample_n(al_dat, size = nrow(al_dat), replace = TRUE)
  # linear model using generalised least squares. The errors are allowed to be correlated and/or have unequal variances.
  am <- am
  #store model coefs
  al_ests[[i]] <- as.list(c(
    smod = coef(am),
    smod.sd.p0 = summary(am)$sigma,
    smod.sd.p = as.numeric(am$modelStruct$varStruct)
  ))
} #end loop

mpars_bs <- bind_cols(
  bind_rows(gr_ests),
  bind_rows(in_ests),
  bind_rows(al_ests)
)

# take the mean of all coefs
mpars_all <- mpars_bs |> 
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

