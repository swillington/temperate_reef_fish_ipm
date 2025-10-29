
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

# Otoliths only using mean parameter values

### run projection functions ----------------------------------------------------------------------------
proj_mean <-
  projectO(
    iparsO = c(l = lower, u = upper, n = n_breaks), #set upper and lower bounds of otolith size
    mpars = apply(mpars_bs, 2, mean),
    ref_temp = mean(unique(useData$covariate)),
    max_a = max_age
  )

## repeat the mean calculations for different temperatures (+/- 1 degree)
proj_mean_m1 <-
  projectO(
    iparsO = c(l = lower, u = upper, n = n_breaks),
    mpars = apply(mpars_bs, 2, mean),
    ref_temp = mean(unique(useData$covariate)) - 1,
    max_a = max_age
  )
proj_mean_p1 <-
  projectO(
    iparsO = c(l = lower, u = upper, n = n_breaks),
    mpars = apply(mpars_bs, 2, mean),
    ref_temp = mean(unique(useData$covariate)) + 1,
    max_a = max_age
  )

### plot the expected age specific density functions -------------------------------------------------------

# create tibble for plotting
proj_mean_df <- as_tibble(proj_mean$n) |> 
  mutate(mesh = proj_mean$mesh) |> 
  select(mesh, starts_with("age_")) |> 
  pivot_longer(cols = starts_with("age_"), 
               names_to = "age", 
               values_to = "density") |>
  mutate(age_num = parse_number(age),     
         age = factor(age, levels = unique(age[order(age_num)])),
         scenario = "mean temp")

# age density plot
ggplot(proj_mean_p1_df, aes(x = mesh, y = density, group = age)) +
  geom_line() +
  coord_cartesian(ylim = c(0, 15)) +
  labs(x = "Otolith Size",
    y = "Density",
    colour = "Age") +
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12))

### plot the otolith size vs age relationship --------------------------------------------------------------------------

# function to calc mean otolith size
summarise_proj_mean <- function(proj_mean_obj, get_quantile_fn) {
  n_mat <- proj_mean_obj$n
  mesh  <- proj_mean_obj$mesh
  dlta  <- proj_mean_obj$dlta
  
  mean_sz      <- apply(n_mat, 2, \(pdens) sum(pdens * mesh * dlta))
  quantile_mat <- apply(n_mat, 2, get_quantile_fn, dlta, mesh)
  
  tibble(
    age     = seq_len(ncol(n_mat)) - 1,
    mean_sz = mean_sz,
    q_low   = quantile_mat[1, ],
    q_high  = quantile_mat[2, ]
  )
}

# create tibble of mean otolith size for plotting
proj_mean_df2 <- summarise_proj_mean(proj_mean, get_quantile) |> 
  mutate(scenario = "Mean temperature")
proj_mean_dfm <- summarise_proj_mean(proj_mean_m1, get_quantile)|> 
  mutate(scenario = "- 1 degree")
proj_mean_dfp <- summarise_proj_mean(proj_mean_p1, get_quantile)|> 
  mutate(scenario = "+ 1 degree")

offset_map <- c("- 1 degree" = -0.2,
                "+ 1 degree" =  0.2,
                "Mean temperature" = 0)

proj_mean_df_all <- bind_rows(proj_mean_df2, proj_mean_dfm, proj_mean_dfp) |> 
  mutate(age_jit = age + offset_map[scenario],
         scenario = factor(scenario, levels = c("Mean temperature", "+ 1 degree", "- 1 degree" ), ordered = TRUE))

# plot expected mean otolith size with age
ggplot(proj_mean_df_all, aes(x = age_jit, y = mean_sz, colour = scenario)) +
  geom_segment(aes(x    = age_jit,
                   xend = age_jit,
                   y    = q_low,
                   yend = q_high,
                   colour = scenario),
              lineend   = "round",
              linewidth = 2, 
              alpha = 0.3) +
  geom_line(linewidth = 0.6) +
  # scale_x_continuous(breaks = unique(summary_df$age),
  #                    labels = unique(summary_df$age)) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("Mean temperature" = "#5fad56", 
                                 "- 1 degree" = "#246eb9", 
                                 "+ 1 degree" = "#f60019"))+
  coord_cartesian(ylim = c(0, 0.7)) +
  labs(
    x = "Age",
    y = "Otolith Size", 
    colour = "Scenario") +
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12))
  
