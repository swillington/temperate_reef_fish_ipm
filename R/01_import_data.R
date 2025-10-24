library(tidyverse)

file_list <- list.files("data/data_raw", full.names = TRUE)

# all data
useData <- read_csv(file_list) %>%
  mutate(covariate = sst_summer) %>% # set temperature that you want to use here
  rename(fishlength = length)


# growth data
gr_dat <-
  useData %>%
  select(
    FishID,
    birth,
    monthcap,
    yearcap,
    annuli_num,
    OL.t1,
    OL.t,
    age,
    covariate)

# allometry data
al_dat <- useData %>% group_by(FishID) %>% 
  slice_max(OL.t, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(rad4 = OL.t) %>% 
  select(c(FishID, birth, fishlength, age, rad4, yearcap, monthcap))
