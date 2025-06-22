formula_cum <- ordinal_value ~ treatment + (1 + treatment | trial)

omr_cum <- brm(
  formula = formula_cum,
  prior = priors,
  sample_prior = TRUE,
  data = all_dat_long,
  family = cumulative(link = "logit"),
  chains = 4,
  cores = 4,
  threads = threading(4),
  backend = "cmdstanr",
  silent = 0,
  seed = seed
)

omr_cum <- add_criterion(omr_cum, "loo")

# Extract the trial-specific coefficients
trial_coefs <- coef(omr_cum, summary = FALSE)$trial[, , "treatmentthrombectomy"] %>%
  posterior_summary() %>%
  as_tibble(rownames = "trial") %>%
  # Sort by Estimate value
  arrange(Estimate) %>%
  # Convert trial to a factor with levels in the sorted order
  mutate(trial = factor(trial, levels = trial))

# Extract the population-level coefficient - corrected approach
pop_coef <- fixef(omr_cum) %>%
  as_tibble(rownames = "Parameter") %>%
  filter(Parameter == "treatmentthrombectomy") %>%
  mutate(trial = "Overall")

# Combine the trial-specific and population-level coefficients
all_coefs <- bind_rows(trial_coefs, pop_coef) %>%
  mutate(trial = factor(trial, levels = c(levels(trial_coefs$trial), "Overall")))

# Create the plot
ggplot(all_coefs, aes(trial, Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_pointrange() +
  coord_flip() +
  labs(x = "trial", y = "Common Odds Ratio") +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  # Add visual distinction for the Overall estimate
  geom_rect(data = subset(all_coefs, trial == "Overall"),
            aes(xmin = as.numeric(trial) - 0.5, 
                xmax = as.numeric(trial) + 0.5,
                ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  # Optional: make the Overall estimate point larger and a different color
  geom_pointrange(data = subset(all_coefs, trial == "Overall"),
                  size = 1.2, color = "red")


sd_parameters <- c("sd_trial__Intercept", "sd_trial__treatmentthrombectomy")
hyp <- paste(sd_parameters, "= 0")
hyp <- hypothesis(omr_cum, hyp, class = NULL)
plot(hyp)

hyp <- hypothesis(omr_cum, "treatmentthrombectomy = 0")
plot(hyp)

hypothesis(
  omr_cum, 
  "cor_trial__Intercept__treatmentthrombectomy = 0",
  class = NULL
) %>% plot()



######

# Extract the trial-specific coefficients
trial_coefs <- coef(omr_cum, summary = FALSE)$trial[, , "treatmentthrombectomy"] %>%
  posterior_summary() %>%
  as_tibble(rownames = "trial") %>%
  # Sort by Estimate value
  arrange(Estimate) %>%
  # Convert trial to a factor with levels in the sorted order
  mutate(trial = factor(trial, levels = trial))

# Extract the population-level coefficient 
pop_coef <- fixef(omr_cum) %>%
  as_tibble(rownames = "Parameter") %>%
  filter(Parameter == "treatmentthrombectomy") %>%
  mutate(trial = "Overall")

# Get posterior draws for predictive interval
posterior_draws <- posterior_samples(omr_cum)
treatment_col <- grep("^b_treatmentthrombectomy$", names(posterior_draws), value = TRUE)
sigma_col <- grep("^sd_trial__treatmentthrombectomy$", names(posterior_draws), value = TRUE)

# Calculate posterior predictive interval
set.seed(123) # For reproducibility
n_draws <- nrow(posterior_draws)
pred_draws <- posterior_draws[, treatment_col] + 
  rnorm(n_draws, 0, posterior_draws[, sigma_col])

# Calculate summary statistics for the predictive interval
pred_summary <- tibble(
  Estimate = mean(pred_draws),
  Q2.5 = quantile(pred_draws, 0.025),
  Q97.5 = quantile(pred_draws, 0.975),
  trial = "Predictive"
)

# Combine all coefficients
all_coefs <- bind_rows(trial_coefs, pop_coef, pred_summary) %>%
  mutate(trial = factor(
    trial, levels = c(levels(trial_coefs$trial), "Overall", "Predictive"))
    ) %>% 
  select(-Parameter)

trial_years <- c(
  "angel" = 2023,
  "attention" = 2022,
  "baoche" = 2022,
  "basics" = 2021,
  "best" = 2020,
  "dawn" = 2018,
  "defuse3" = 2018,
  "distal" = 2025,
  "escape" = 2015,
  "escapemevo" = 2025,
  "extend" = 2015,
  "ims3" = 2013,
  "mrclean" = 2015,
  "mrcleanlate" = 2023,
  "mrrescue" = 2013,
  "piste" = 2017, 
  "positive" = 2022,
  "rescue" = 2022,
  "resilient" = 2020,
  "revascat" = 2015,
  "select" = 2023,
  "swift" = 2015, 
  "synthesis" = 2013,
  "tension" = 2023,
  "tesla" = 2024,
  "therapy" = 2016,
  "thrace" = 2016,
  "Overall" = 3000,
  "Predictive" = 3001
)

# 1) Assign a 'year' column
all_coefs <- all_coefs %>%
  mutate(year = trial_years[as.character(trial)])  

# 2) Sort levels of 'trial' by year; earliest year = top row after coord_flip
#    If you want the earliest at the *bottom*, swap out decreasing=FALSE -> decreasing=TRUE
all_coefs <- all_coefs %>%
  mutate(trial = factor(trial,
                        levels = names(sort(trial_years, decreasing = TRUE))))

# 3) Re-plot with newly factored 'trial'
ggplot(all_coefs, aes(trial, Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_pointrange() +
  coord_flip() +
  labs(x = "Trial", y = "Common Odds Ratio") +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  # Add the highlight rectangle for Overall/Predictive
  geom_rect(data = subset(all_coefs, trial %in% c("Overall", "Predictive")),
            aes(xmin = as.numeric(trial) - 0.5, 
                xmax = as.numeric(trial) + 0.5,
                ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3) +
  # Larger, red point for Overall
  geom_pointrange(data = subset(all_coefs, trial == "Overall"),
                  size = 0.5, color = "red") +
  # Larger, blue (dashed) point for Predictive
  geom_pointrange(data = subset(all_coefs, trial == "Predictive"),
                  size = 0.5, color = "blue", linetype = "dashed")

#### posterior difference

new_data <- data.frame(treatment = c("medical", "thrombectomy"))
post_cat <- posterior_epred(omr_cum, newdata = new_data, re_formula = NA)
diff_mat <- post_cat[,2,] - post_cat[,1,]  # S x C matrix
diff_summary <- apply(diff_mat, 2, 
                      function(x) c(mean(x),
                                    quantile(x, 0.025),
                                    quantile(x, 0.975)))
df_diff <- data.frame(
  # Change here: subtract 1 from each category value to shift from 1:7 to 0:6
  category   = factor(seq_len(ncol(diff_mat)) - 1),
  diff_mean  = diff_summary[1, ],
  diff_lower = diff_summary[2, ],
  diff_upper = diff_summary[3, ]
)
ggplot(df_diff, aes(x = category, y = diff_mean)) +
  geom_point() +
  geom_linerange(aes(ymin = diff_lower, ymax = diff_upper)) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(
    x = "mRS Category", 
    y = "Posterior mean difference (Thrombectomy - Medical)"
  ) +
  theme_bw()

#### ordinal forest plot

# 1) Build a data frame for each (Trial x Treatment) combination
#    This ensures we get separate predictions for medical vs thrombectomy 
#    *in each trial*, incorporating that trial's random effect.
df_trial <- expand.grid(
  trial     = unique(all_dat_long$trial),
  treatment = c("medical", "thrombectomy")
)
# If your model has other covariates, you must include them here too 
# (e.g., center them or pick reference values, etc.)

# 2) Posterior predictions WITH random effects
#    re_formula = NULL => Use the random effects by trial
post_all <- posterior_epred(omr_cum, newdata = df_trial, re_formula = NULL)
# Dimensions of returned array: S x N x C
#    S = posterior draws
#    N = number of rows in df_trial
#    C = number of ordinal categories

# Add a row index so we can keep track
df_trial <- df_trial %>% mutate(row_index = row_number())

# 3) Pair up medical vs thrombectomy within each trial and compute differences
df_pairs <- df_trial %>%
  group_by(trial) %>%
  summarize(
    med_idx = row_index[treatment == "medical"],
    thr_idx = row_index[treatment == "thrombectomy"]
  ) %>%
  ungroup()

# Loop over trials, compute posterior difference for each ordinal category
all_diffs <- vector("list", nrow(df_pairs))

for (i in seq_len(nrow(df_pairs))) {
  
  idx_med <- df_pairs$med_idx[i]
  idx_thr <- df_pairs$thr_idx[i]
  
  # diff_mat is S x C
  diff_mat <- post_all[, idx_thr, ] - post_all[, idx_med, ]
  
  # Summaries across posterior draws
  diff_summary <- apply(diff_mat, 2, function(x) {
    c(mean = mean(x),
      lower = quantile(x, 0.025),
      upper = quantile(x, 0.975))
  })
  
  df_i <- data.frame(
    trial = df_pairs$trial[i],
    # Shift category from 1..7 to 0..6 for typical mRS labeling
    category = factor(seq_len(ncol(diff_mat)) - 1),
    diff_mean  = diff_summary[1, ],
    diff_lower = diff_summary[2, ],
    diff_upper = diff_summary[3, ]
  )
  
  all_diffs[[i]] <- df_i
}

df_final <- do.call(rbind, all_diffs)

# 4) Do the same for "overall" (population-level) effect
#    We'll set re_formula=NA to EXCLUDE random effects => population-level
df_overall <- data.frame(treatment = c("medical", "thrombectomy"))
post_overall <- posterior_epred(omr_cum, newdata = df_overall, re_formula = NA)

diff_mat_overall <- post_overall[, 2, ] - post_overall[, 1, ]
diff_summary_overall <- apply(diff_mat_overall, 2, function(x) {
  c(mean = mean(x),
    lower = quantile(x, 0.025),
    upper = quantile(x, 0.975))
})
df_overall_final <- data.frame(
  trial = "Overall",
  category = factor(seq_len(ncol(diff_mat_overall)) - 1),
  diff_mean  = diff_summary_overall[1, ],
  diff_lower = diff_summary_overall[2, ],
  diff_upper = diff_summary_overall[3, ]
)

# Combine
df_plot <- bind_rows(df_final, df_overall_final)

trial_years <- c(
  "angel" = 2023,
  "attention" = 2022,
  "baoche" = 2022,
  "basics" = 2021,
  "best" = 2020,
  "dawn" = 2018,
  "defuse3" = 2018,
  "distal" = 2025,
  "escape" = 2015,
  "escapemevo" = 2025,
  "extend" = 2015,
  "ims3" = 2013,
  "mrclean" = 2015,
  "mrcleanlate" = 2023,
  "mrrescue" = 2013,
  "piste" = 2017, 
  "positive" = 2022,
  "rescue" = 2022,
  "resilient" = 2020,
  "revascat" = 2015,
  "select" = 2023,
  "swift" = 2015, 
  "synthesis" = 2013,
  "tension" = 2023,
  "tesla" = 2024,
  "therapy" = 2016,
  "thrace" = 2016,
  "Overall" = Inf
)

# 1. Add a 'year' column
df_plot <- df_plot %>%
  mutate(year = trial_years[trial])

# 2. Reorder 'trial' factor by ascending 'year'
df_plot <- df_plot %>%
  mutate(trial = reorder(trial, year))

# Now facet_wrap will follow the trial factor levels (i.e. by year)
ggplot(df_plot, aes(x = category, y = diff_mean)) +
  geom_point(size = 0.5) +
  geom_linerange(aes(ymin = diff_lower, ymax = diff_upper)) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~ trial, ncol = 7) +
  labs(
    x = "mRS Category",
    y = "Posterior Mean Difference (Thrombectomy − Medical)"
  ) +
  theme_bw()

