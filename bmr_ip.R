library(dplyr)

# Bernoulli Meta-Regression with Informed Priors
all_dat <- read_rds(here("all_dat.rds"))

dichot_dat <- all_dat %>%
  rowwise() %>%
  mutate(
    favorable = sum(c(zero, one, two), na.rm = TRUE),
    unfavorable = sum(c(three, four, five, six), na.rm = TRUE)
  ) %>%
  # Select only the columns we want to keep
  select(type, type_id, trial, trial_id, treatment, treatment_id, 
         favorable, unfavorable)
  
print(dichot_dat, n = Inf)

dichot_dat_long <- dichot_dat %>%
  pivot_longer(
    cols = -c(type, type_id, trial, trial_id, treatment, treatment_id),
    names_to = "ordinal_value", 
    values_to = "count"
  ) %>% 
  uncount(count) %>%
  arrange(treatment_id, ordinal_value) %>% 
  mutate(ordinal_value = case_when(
    ordinal_value == "favorable" ~ 1,
    ordinal_value == "unfavorable" ~ 0
  ))

formula <- ordinal_value ~ 0 + type + treatment + (1 + treatment | trial)
priors <- c(prior(normal(0, 1), class = "b"),
            prior(normal(0, 1), class = "sd"))

bmr_ip <- brm(
  formula = formula,
  prior = priors, 
  data = dichot_dat_long,
  family = bernoulli(link = "logit"),
  chains = 4, 
  cores = 4, 
  backend = "cmdstanr",
  # refresh = 0,
  seed = seed
)

bmr_ip <- add_criterion(bmr_ip, "loo")

saveRDS(bmr_ip, here("bmr_ip.rds"))

# calculate posterior absolute risk contrast by stroke type

bmr_absolute_contrast <- function(model, types = c("large", "early", "late", "basilar")) {
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
    diff_posterior <- pp_treat - pp_control
    
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

bmr_results <- bmr_absolute_contrast(bmr_ip) %>% 
  relocate(type, .before = mean)

print(bmr_results)
