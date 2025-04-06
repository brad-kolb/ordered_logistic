# load packages
library(tidyverse)
library(brms)
library(cmdstanr)
library(here)

#### import and prepare data ####

seed <- 123

all_dat <- read.csv(here("thrombectomy.csv")) %>% 
  as_tibble() %>% 
  rename(type = group) %>% 
  rename(type_id = group_id) %>% 
  mutate(treatment = recode(treatment,
                            "treatment" = "thrombectomy", 
                            "control" = "medical")) %>% 
  # tension trial collapses mrs 0 and 1
  # "mRS  scores  0  and  1  were  merged  due  to  the  very  low   numbers  
  # in  the  best  medical treatment  group."
  # 10 in treatment group were 0 or 1
  # 2 in medical group were 0 or 1
  # we divided them evenly between the two groups
  mutate(one = replace_na(one, 0))

# change to factors
all_dat$type <- factor(all_dat$type, 
                       levels = unique(all_dat$type[order(all_dat$type_id)]))
all_dat$trial <- factor(all_dat$trial, 
                        levels = unique(all_dat$trial[order(all_dat$trial_id)]))
all_dat$treatment <- factor(all_dat$treatment,
                            levels = unique(all_dat$treatment[order(all_dat$treatment_id)]))

print(all_dat, n = Inf)

saveRDS(all_dat, here("all_dat.rds"))

#### fit model ####
all_dat_long <- all_dat %>%
  pivot_longer(
    cols = -c(type, type_id, trial, trial_id, treatment, treatment_id),
    names_to = "ordinal_value", 
    values_to = "count"
  ) %>% 
  uncount(count) %>%
  arrange(treatment_id, ordinal_value) %>% 
  mutate(ordinal_value = case_when(
    ordinal_value == "zero" ~ 1,
    ordinal_value == "one" ~ 2, 
    ordinal_value == "two" ~ 3,
    ordinal_value == "three" ~ 4,
    ordinal_value == "four" ~ 5,
    ordinal_value == "five" ~ 6,
    ordinal_value == "six" ~ 7
  ))

formula <- ordinal_value ~ type + cs(treatment) + (1 + treatment | trial)
priors <- c(prior(normal(0, 1), class = "b"),
            prior(normal(0, 1), class = "sd"),
            prior(normal(0, 1), class = "Intercept"))


omr_ip <- brm(
  formula = formula,
  prior = priors, 
  data = all_dat_long,
  family = acat(link = "logit"),
  chains = 4, 
  cores = 4, 
  backend = "cmdstanr",
  # refresh = 0,
  seed = seed
)

omr_ip <- add_criterion(omr_ip, "loo")

saveRDS(omr_ip, here("omr_ip.rds"))

# absolute contrasts
categories <- c("large", "early", "late", "basilar", "medium", "historical")

omr_ip_ac <- function(model, types = categories) {
  # Initialize an empty list to store results
  results <- list()
  
  # Loop through each type
  for (current_type in types) {
    # Create new data frames for each type
    newdata_control <- data.frame(treatment = "medical", type = current_type)
    newdata_treat <- data.frame(treatment = "thrombectomy", type = current_type)
    
    # Calculate posterior predictions
    pp_control <- posterior_epred(
      model,
      newdata = newdata_control,
      re_formula = NA,
      ndraws = NULL
    )
    
    pp_treat <- posterior_epred(
      model,
      newdata = newdata_treat,
      re_formula = NA
    )
    
    # Calculate difference
    diff_posterior <- rowSums(pp_treat[, 1, 1:3]) - rowSums(pp_control[, 1, 1:3])
    
    # Store summary for current type
    results[[current_type]] <- tibble(
      mean = mean(diff_posterior),
      "2.5%" = quantile(diff_posterior, 0.025),
      "97.5%" = quantile(diff_posterior, 0.975)
    ) %>%
      mutate(type = current_type)
  }
  
  # Combine all results into a single data frame
  final_summary <- bind_rows(results) 
  
  return(final_summary)
}

omr_results <- omr_ip_ac(omr_ip) %>% 
  relocate(type, .before = mean)

print(omr_results)

#### sensitivity checks ####

# flat priors
omr_flat <- brm(
  formula = formula,
  data = all_dat_long,
  family = acat(link = "logit"),
  chains = 4,
  cores = 4,
  backend = "cmdstanr",
  # refresh = 0,
  seed = seed
)

omr_flat <- add_criterion(omr_flat, "loo")

saveRDS(omr_flat, here("omr_flat.rds"))

# proportional odds model

formula_cum <- ordinal_value ~ type + treatment + (1 + treatment | trial)

omr_cum <- brm(
  formula = formula_cum,
  prior = priors,
  data = all_dat_long,
  family = cumulative(link = "logit"),
  chains = 4,
  cores = 4,
  backend = "cmdstanr",
  # refresh = 0,
  seed = seed
)

omr_cum <- add_criterion(omr_cum, "loo")

saveRDS(omr_cum, here("omr_cum.rds"))

# drop stroke type predictors

formula_oma <- ordinal_value ~ cs(treatment) + (1 + treatment | trial)

oma <- brm(
  formula = formula_oma,
  prior = priors,
  data = all_dat_long,
  family = acat(link = "logit"),
  chains = 4,
  cores = 4,
  backend = "cmdstanr",
  # refresh = 0,
  seed = seed
)

oma <- add_criterion(oma, "loo")

saveRDS(oma, here("oma.rds"))
