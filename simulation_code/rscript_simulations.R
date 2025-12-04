#!/usr/bin/env Rscript

# Command-line interface for run_full_simulation_study
# Usage: Rscript run_simulation_cli.R --source_file ./simulation_functions.R --inner_interval 0.90 --n_replicates 50 --output_dir ./results --output_name simulation_results.rds

library(optparse)
library(orthomr)
library(polynom)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(TwoSampleMR)

# Define command-line options
option_list <- list(
  make_option(c("--source_file"), type="character", default=NULL,
              help="Path to R file containing simulation functions (required)"),
  make_option(c("--inner_interval"), type="double", default=0.90,
              help="Proportion of exposure distribution to use for error calculation [default %default]"),
  make_option(c("--n_replicates"), type="integer", default=50,
              help="Number of simulations per scenario [default %default]"),
  make_option(c("--prop_int"), type="character", default="0,0.1,0.25,0.5",
              help="Comma-separated proportion of pleiotropic variants [default %default]"),
  make_option(c("--inst_int"), type="character", default="30,50",
              help="Comma-separated number of instruments [default %default]"),
  make_option(c("--ss_int"), type="character", default="10000,25000,50000",
              help="Comma-separated sample sizes [default %default]"),
  make_option(c("--output_dir"), type="character", default=".",
              help="Directory to save output RDS file [default %default]"),
  make_option(c("--output_name"), type="character", default="simulation_results.rds",
              help="Name of output RDS file [default %default]")
)

# Parse arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#' Run Usual MR
#'
#' @param sim_data Simulated data from sim_exp_out_non_linear
#'
#' @return Data frame with IVW MR results
run_ivw <- function(sim_data){
  gwas_temp <- data.frame()
  for (i in 1:ncol(sim_data$genotypes)) {
    model <- lm(sim_data$simulated_data$exposure_observed ~ sim_data$genotypes[, i])
    gwas_temp <- rbind(gwas_temp, data.frame(
      ins = paste0("X", i),
      beta = coefficients(model)[2],
      se = coef(summary(model))[, "Std. Error"][2],
      pval = coef(summary(model))[, "Pr(>|t|)"][2],
      af = sim_data$allele_frequencies[i]
    ))
  }
  
  gwas_filt<-gwas_temp %>%
    filter(pval<5e-8) %>%
    mutate(EA="C") %>%
    mutate(OA="A")
  
  gwas_temp_out <- data.frame()
  for (i in 1:ncol(sim_data$genotypes)) {
    model <- lm(sim_data$simulated_data$outcome_observed ~ sim_data$genotypes[, i])
    gwas_temp_out <- rbind(gwas_temp_out, data.frame(
      ins = paste0("X", i),
      beta = coefficients(model)[2],
      se = coef(summary(model))[, "Std. Error"][2],
      pval = coef(summary(model))[, "Pr(>|t|)"][2],
      af = sim_data$allele_frequencies[i]
    ))
  }
  
  gwas_filt_out<-gwas_temp_out %>%
    filter(ins %in% gwas_filt$ins) %>%
    mutate(EA="C") %>%
    mutate(OA="A")
  
  mr_exp<-format_data(dat=gwas_filt,
                      type="exposure",
                      snp_col="ins",
                      beta_col="beta",
                      se_col="se",
                      pval_col="pval",
                      eaf_col="af",
                      effect_allele_col = "EA",
                      other_allele_col = "OA")
  
  mr_out<-format_data(dat=gwas_filt_out,
                      type="outcome",
                      snp_col="ins",
                      beta_col="beta",
                      se_col="se",
                      pval_col="pval",
                      eaf_col="af",
                      effect_allele_col = "EA",
                      other_allele_col = "OA")
  
  mr_harm<-harmonise_data(mr_exp,mr_out)
  
  mr_res<-mr(mr_harm,method_list=c("mr_ivw"))
  
  return(mr_res)
}


#' Run Single Simulation Scenario
#'
#' @param scenario_name Name of scenario
#' @param true_func True causal function
#' @param sample_size Sample size
#' @param n_instruments Number of genetic instruments
#' @param max_degree Maximum polynomial degree to fit
#' @param include_pleiotropy Include pleiotropic effects
#' @param prop_pleiotropic Proportion of pleiotropic variants
#' @param prop_variance_explained Proportion of outcome variance explained
#' @param n_replicates Number of simulation replicates
#' @param inner_interval Proportion of exposure distribution to use for error calculation (default: 0.90)
#' @param compare_ivw If TRUE, also run standard IVW MR for comparison (default: TRUE)
#'
#' @return Data frame with performance metrics
run_simulation_scenario <- function(scenario_name,
                                    true_func,
                                    sample_size = 10000,
                                    n_instruments = 30,
                                    max_degree = 5,
                                    include_pleiotropy = FALSE,
                                    prop_pleiotropic = 0.25,
                                    prop_variance_explained = 0.5,
                                    n_replicates = 20,
                                    inner_interval = 0.90,
                                    compare_ivw = TRUE) {
  
  results_list <- list()
  
  for (rep in 1:n_replicates) {
    
    # Set seed for reproducibility
    #set.seed(1000 + rep)
    
    # Simulate data
    sim_data <- sim_exp_out_non_linear(
      sample_size = sample_size,
      n_instruments = n_instruments,
      exposure_to_outcome_func = true_func,
      include_pleiotropy = include_pleiotropy,
      prop_pleiotropic = prop_pleiotropic,
      prop_variance_explained = prop_variance_explained
    )
    
    # Extract variables
    outcome_observed <- sim_data[["simulated_data"]]$outcome_observed
    exposure_observed <- sim_data[["simulated_data"]]$exposure_observed
    outcome_true <- sim_data[["simulated_data"]]$outcome_true
    exposure_true <- sim_data[["simulated_data"]]$exposure_true
    genotypes <- sim_data[["genotypes"]]
    allele_frequencies <- sim_data[["allele_frequencies"]]
    
    # Create orthogonal polynomials
    poly_result <- orthopol(degree = max_degree, exposure_values = exposure_observed)
    
    # Run GWAS for exposures
    gwas_exposure_list <- list()
    for (j in 1:max_degree) {
      gwas_temp <- data.frame()
      for (i in 1:ncol(genotypes)) {
        model <- lm(poly_result[["values"]][, j] ~ genotypes[, i])
        gwas_temp <- rbind(gwas_temp, data.frame(
          ins = paste0("X", i),
          beta = coefficients(model)[2],
          se = coef(summary(model))[, "Std. Error"][2],
          pval = coef(summary(model))[, "Pr(>|t|)"][2],
          af = allele_frequencies[i]
        ))
      }
      gwas_exposure_list[[j]] <- gwas_temp
    }
    
    # Run GWAS for outcome
    gwas_outcome <- data.frame()
    for (i in 1:ncol(genotypes)) {
      model <- lm(outcome_observed ~ genotypes[, i])
      gwas_outcome <- rbind(gwas_outcome, data.frame(
        ins = paste0("X", i),
        beta = coefficients(model)[2],
        se = coef(summary(model))[, "Std. Error"][2],
        pval = coef(summary(model))[, "Pr(>|t|)"][2],
        af = allele_frequencies[i]
      ))
    }
    
    # Format for MVMR
    mvmr_exposure_data <- merge_exp_inst(gwas_exposure_list)
    
    mvmr_outcome_data <- gwas_outcome %>%
      filter(ins %in% mvmr_exposure_data$SNP) %>%
      rename(SNP = ins) %>%
      mutate(
        outcome = "outcome",
        id.outcome = "outcome",
        effect_allele.outcome = "A",
        other_allele.outcome = "C"
      ) %>%
      rename(
        beta.outcome = beta,
        se.outcome = se,
        pval.outcome = pval,
        eaf.outcome = af
      )
    
    # Run MR analysis (with error handling)
    tryCatch({
      analysis_results <- mr_sim_res(
        mvmr_exposure_data = mvmr_exposure_data,
        mvmr_outcome_data = mvmr_outcome_data,
        simulation_data = sim_data,
        inner_interval = 0.90,
        include_intercept = TRUE,
        orthogonal_poly_coef = poly_result$coefficients
      )
      
      # Calculate metrics
      # Limit to inner interval of exposure distribution
      lower_quantile <- (1 - inner_interval) / 2
      upper_quantile <- 1 - lower_quantile
      quantiles <- quantile(exposure_observed, c(lower_quantile, upper_quantile))
      
      analysis_df <- data.frame(
        exposure_true = exposure_true,
        exposure_observed = exposure_observed,
        outcome_true = outcome_true,
        outcome_observed = outcome_observed
      ) %>%
        filter(exposure_observed >= quantiles[1] & exposure_observed <= quantiles[2])
      
      # Also calculate metrics on full data for comparison
      full_df <- data.frame(
        exposure_true = exposure_true,
        exposure_observed = exposure_observed,
        outcome_true = outcome_true,
        outcome_observed = outcome_observed
      )
      
      # Predict outcomes (inner interval)
      predicted_outcome <- predict_mvmr(
        analysis_df$exposure_true,
        analysis_results$polynomial_coefficients
      )
      
      # Predict outcomes (full data)
      predicted_outcome_full <- predict_mvmr(
        full_df$exposure_true,
        analysis_results$polynomial_coefficients
      )
      
      # Calculate metrics on inner interval
      mse <- mean((analysis_df$outcome_true - predicted_outcome)^2)
      mae <- mean(abs(analysis_df$outcome_true - predicted_outcome))
      bias <- mean(analysis_df$outcome_true - predicted_outcome)
      r_squared <- cor(analysis_df$outcome_true, predicted_outcome)^2
      
      # Calculate metrics on full data
      mse_full <- mean((full_df$outcome_true - predicted_outcome_full)^2)
      mae_full <- mean(abs(full_df$outcome_true - predicted_outcome_full))
      bias_full <- mean(full_df$outcome_true - predicted_outcome_full)
      r_squared_full <- cor(full_df$outcome_true, predicted_outcome_full)^2
      
      # IVW predictions (linear only)
      ivw_mse <- NA
      ivw_mae <- NA
      ivw_bias <- NA
      ivw_r_squared <- NA
      ivw_mse_full <- NA
      ivw_mae_full <- NA
      ivw_bias_full <- NA
      ivw_r_squared_full <- NA
      ivw_beta <- NA
      ivw_se <- NA
      ivw_pval <- NA
      
      if (compare_ivw) {
        ivw_results<-run_ivw(sim_data)
        
        ivw_beta <- ivw_results$b[1]
        ivw_se <- ivw_results$se[1]
        ivw_pval <- ivw_results$pval[1]
        
        # Predict using linear IVW estimate
        predicted_ivw <- ivw_beta * analysis_df$exposure_true
        predicted_ivw_full <- ivw_beta * full_df$exposure_true
        
        # Calculate IVW metrics (inner interval)
        ivw_mse <- mean((analysis_df$outcome_true - predicted_ivw)^2)
        ivw_mae <- mean(abs(analysis_df$outcome_true - predicted_ivw))
        ivw_bias <- mean(analysis_df$outcome_true - predicted_ivw)
        ivw_r_squared <- cor(analysis_df$outcome_true, predicted_ivw)^2
        
        # Calculate IVW metrics (full data)
        ivw_mse_full <- mean((full_df$outcome_true - predicted_ivw_full)^2)
        ivw_mae_full <- mean(abs(full_df$outcome_true - predicted_ivw_full))
        ivw_bias_full <- mean(full_df$outcome_true - predicted_ivw_full)
        ivw_r_squared_full <- cor(full_df$outcome_true, predicted_ivw_full)^2
        
        
        
      }
      
      # Calculate bias at different exposure quantiles
      exposure_quintiles <- quantile(analysis_df$exposure_observed, seq(0, 1, 0.2))
      bias_by_quintile <- sapply(1:(length(exposure_quintiles)-1), function(q) {
        idx <- which(analysis_df$exposure_observed >= exposure_quintiles[q] &
                       analysis_df$exposure_observed < exposure_quintiles[q+1])
        if(length(idx) > 0) {
          mean(analysis_df$outcome_true[idx] - predicted_outcome[idx])
        } else {
          NA
        }
      })
      
      # Determine selected degree (count non-zero significant terms)
      selected_degree <- nrow(analysis_results$mvmr_results)
      if("Intercept" %in% analysis_results$mvmr_results$exposure) {
        selected_degree <- selected_degree - 1
      }
      
      # Store results
      results_list[[rep]] <- data.frame(
        scenario = scenario_name,
        replicate = rep,
        sample_size = sample_size,
        n_instruments = n_instruments,
        inner_interval = inner_interval,
        # Polynomial MVMR: Inner interval metrics
        poly_mse = mse,
        poly_mae = mae,
        poly_bias = bias,
        poly_r_squared = r_squared,
        # Polynomial MVMR: Full data metrics
        poly_mse_full = mse_full,
        poly_mae_full = mae_full,
        poly_bias_full = bias_full,
        poly_r_squared_full = r_squared_full,
        # IVW: Inner interval metrics
        ivw_mse = ivw_mse,
        ivw_mae = ivw_mae,
        ivw_bias = ivw_bias,
        ivw_r_squared = ivw_r_squared,
        # IVW: Full data metrics
        ivw_mse_full = ivw_mse_full,
        ivw_mae_full = ivw_mae_full,
        ivw_bias_full = ivw_bias_full,
        ivw_r_squared_full = ivw_r_squared_full,
        # IVW estimates
        ivw_beta = ivw_beta,
        ivw_se = ivw_se,
        ivw_pval = ivw_pval,
        # Other metrics
        selected_degree = selected_degree,
        bias_q1 = bias_by_quintile[1],
        bias_q2 = bias_by_quintile[2],
        bias_q3 = bias_by_quintile[3],
        bias_q4 = bias_by_quintile[4],
        bias_q5 = bias_by_quintile[5],
        min_f_stat = min(analysis_results$instrument_strength),
        n_samples_inner = nrow(analysis_df),
        n_samples_full = nrow(full_df),
        converged = TRUE
      )
      
    }, error = function(e) {
      results_list[[rep]] <- data.frame(
        scenario = scenario_name,
        replicate = rep,
        sample_size = sample_size,
        n_instruments = n_instruments,
        inner_interval = inner_interval,
        poly_mse = NA,
        poly_mae = NA,
        poly_bias = NA,
        poly_r_squared = NA,
        poly_mse_full = NA,
        poly_mae_full = NA,
        poly_bias_full = NA,
        poly_r_squared_full = NA,
        ivw_mse = NA,
        ivw_mae = NA,
        ivw_bias = NA,
        ivw_r_squared = NA,
        ivw_mse_full = NA,
        ivw_mae_full = NA,
        ivw_bias_full = NA,
        ivw_r_squared_full = NA,
        ivw_beta = NA,
        ivw_se = NA,
        ivw_pval = NA,
        selected_degree = NA,
        bias_q1 = NA,
        bias_q2 = NA,
        bias_q3 = NA,
        bias_q4 = NA,
        bias_q5 = NA,
        min_f_stat = NA,
        n_samples_inner = NA,
        n_samples_full = NA,
        converged = FALSE
      )
    })
    
    if(rep %% 10 == 0) {
      message(paste0("Scenario: ", scenario_name, " - Replicate: ", rep, "/", n_replicates))
    }
  }
  
  # Combine results
  return(bind_rows(results_list))
}

#' Run Full Simulation Study
#'
#' @param inner_interval Proportion of exposure distribution to use for error calculation (default: 0.90)
#' @param n_replicates Number of simulations per scenario (default: 50)
#'
#' @return Data frame with all simulation results
run_full_simulation_study <- function(inner_interval = 0.90,
                                      n_replicates=50,
                                      prop_int=c(0, 0.1, 0.25, 0.5),
                                      inst_int=c(30, 50),
                                      ss_int=c(10000,25000,50000)) {
  
  all_results <- list()
  for(prop in prop_int) {
    for(n_inst in inst_int) {
      for(ss in ss_int){
        message(paste("Running with sample size:", ss))
        message(paste("Running with pleiotropy prop:", prop))
        message(paste("Running with n_instruments:", n_inst))
        for(n in 1:n_replicates){
          
          cust_thresh<-function(x){
            cutoff<-runif(1,-10,2)
            plateau<-runif(1,-10,10)
            slope<-runif(1,-10,10)
            df<-data.frame(y_bool=(x>cutoff)) %>%
              mutate(y_tmp=plateau+(x-cutoff)*slope) %>%
              mutate(y=ifelse(y_bool, y_tmp, plateau))
            return(df %>% pull(y))
          }
          
          scenarios <- list(
            # Functional forms
            linear = list(
              name = "Linear",
              func = as.function(polynomial(coef = c(runif(1,-10,10),
                                                     runif(1,-10,10)))),
              true_degree = 1
            ),
            quadratic = list(
              name = "Quadratic",
              func = as.function(polynomial(coef = c(runif(1,-10,10),
                                                     runif(1,-10,10),
                                                     runif(1,-10,10)))),
              true_degree = 2
            ),
            cubic = list(
              name = "Cubic",
              func = as.function(polynomial(coef = c(runif(1,-10,10),
                                                     runif(1,-10,10),
                                                     runif(1,-10,10),
                                                     runif(1,-10,10)))),
              true_degree = 3
            ),
            fourth_order = list(
              name = "Quartic",
              func = as.function(polynomial(coef = c(runif(1,-10,10),
                                                     runif(1,-10,10),
                                                     runif(1,-10,10),
                                                     runif(1,-10,10),
                                                     runif(1,-10,10)))),
              true_degree = 4
            ),
            threshold = list(
              name = "Threshold",
              func = cust_thresh,
              true_degree = NA
            )
          )
          
          for(scenario in names(scenarios)) {
            result<-c()
            message(paste0("Running scenario: ", scenarios[[scenario]]$name, ", pleiotropy ", prop*100,"%, n_instruments: ",n_inst, ", ss: ",ss, ", replicate: ",n))
            result <- tryCatch({
              run_simulation_scenario(
                scenario_name = paste0(paste0(scenarios[[scenario]]$name,"_pleiotropy_",prop,"_n_ins_",n_inst,"_ss_",ss,"_rep_",n)),
                true_func = scenarios[[scenario]]$func,
                include_pleiotropy = (prop > 0),
                prop_pleiotropic = prop,
                sample_size = ss,
                n_instruments = n_inst,
                n_replicates = 1,
                inner_interval = inner_interval,
                compare_ivw = TRUE
              ) %>%
                mutate(prop_pleio=prop) %>%
                bind_rows(result,.)
            }, error = function(e) {
              message(paste0("ERROR in scenario", scenarios[[scenario]]$name,"_pleiotropy_",prop,"_n_ins_",n_inst, "_ss_", ss, e$message))
              return(NULL)
            })
            
            if (!is.null(result) && nrow(result) > 0) {
              all_results[[paste0(scenarios[[scenario]]$name,"_pleiotropy_",prop,"_n_ins_",n_inst,"_ss_",ss,"_rep_",n)]] <- result
              message(paste("  -> Success:", nrow(result), "rows"))
            }
          }
        }
      }
    }
  }
  
  # Combine all results
  message("=== Combining Results ===")
  if (length(all_results) == 0) {
    warning("No successful simulation scenarios! All failed.")
    return(data.frame())
  }
  
  final_results <- bind_rows(all_results)
  message(paste("Total rows in final results:", nrow(final_results)))
  
  return(final_results)
}

#' Summarize Simulation Results
#'
#' @param results Data frame from run_full_simulation_study()
#' @param use_full_data If TRUE, summarize full data metrics; if FALSE, summarize inner interval metrics (default: FALSE)
#' @param method Method to summarize: "polynomial" or "ivw" (default: "polynomial")
#'
#' @return Summary statistics by scenario
summarize_results <- function(results, use_full_data = FALSE, method = "polynomial") {
  
  # Select appropriate columns based on method and data type
  if (method == "polynomial") {
    if(use_full_data) {
      results %>%
        filter(converged == TRUE) %>%
        mutate(scenario=stringr::str_replace(scenario,"_rep_[0-9]*","")) %>%
        group_by(scenario) %>%
        summarise(
          n_converged = n(),
          mean_mse = mean(poly_mse_full, na.rm = TRUE),
          sd_mse = sd(poly_mse_full, na.rm = TRUE),
          mean_mae = mean(poly_mae_full, na.rm = TRUE),
          mean_bias = mean(poly_bias_full, na.rm = TRUE),
          sd_bias = sd(poly_bias_full, na.rm = TRUE),
          mean_r2 = mean(poly_r_squared_full, na.rm = TRUE),
          median_selected_degree = median(selected_degree, na.rm = TRUE),
          mean_min_f = mean(min_f_stat, na.rm = TRUE),
          method = "Polynomial",
          data_type = "Full",
          .groups = "drop"
        ) %>%
        arrange(mean_mse)
    } else {
      results %>%
        filter(converged == TRUE) %>%
        mutate(scenario=stringr::str_replace(scenario,"_rep_[0-9]*","")) %>%
        group_by(scenario) %>%
        summarise(
          n_converged = n(),
          mean_mse = mean(poly_mse, na.rm = TRUE),
          sd_mse = sd(poly_mse, na.rm = TRUE),
          mean_mae = mean(poly_mae, na.rm = TRUE),
          mean_bias = mean(poly_bias, na.rm = TRUE),
          sd_bias = sd(poly_bias, na.rm = TRUE),
          mean_r2 = mean(poly_r_squared, na.rm = TRUE),
          median_selected_degree = median(selected_degree, na.rm = TRUE),
          mean_min_f = mean(min_f_stat, na.rm = TRUE),
          inner_interval = first(inner_interval),
          method = "Polynomial",
          data_type = "Inner",
          .groups = "drop"
        ) %>%
        arrange(mean_mse)
    }
  } else if (method == "ivw") {
    if(use_full_data) {
      results %>%
        filter(converged == TRUE) %>%
        mutate(scenario=stringr::str_replace(scenario,"_rep_[0-9]*","")) %>%
        group_by(scenario) %>%
        summarise(
          n_converged = n(),
          mean_mse = mean(ivw_mse_full, na.rm = TRUE),
          sd_mse = sd(ivw_mse_full, na.rm = TRUE),
          mean_mae = mean(ivw_mae_full, na.rm = TRUE),
          mean_bias = mean(ivw_bias_full, na.rm = TRUE),
          sd_bias = sd(ivw_bias_full, na.rm = TRUE),
          mean_r2 = mean(ivw_r_squared_full, na.rm = TRUE),
          mean_ivw_beta = mean(ivw_beta, na.rm = TRUE),
          mean_min_f = mean(min_f_stat, na.rm = TRUE),
          method = "IVW",
          data_type = "Full",
          .groups = "drop"
        ) %>%
        arrange(mean_mse)
    } else {
      results %>%
        filter(converged == TRUE) %>%
        mutate(scenario=stringr::str_replace(scenario,"_rep_[0-9]*","")) %>%
        group_by(scenario) %>%
        summarise(
          n_converged = n(),
          mean_mse = mean(ivw_mse, na.rm = TRUE),
          sd_mse = sd(ivw_mse, na.rm = TRUE),
          mean_mae = mean(ivw_mae, na.rm = TRUE),
          mean_bias = mean(ivw_bias, na.rm = TRUE),
          sd_bias = sd(ivw_bias, na.rm = TRUE),
          mean_r2 = mean(ivw_r_squared, na.rm = TRUE),
          mean_ivw_beta = mean(ivw_beta, na.rm = TRUE),
          inner_interval = first(inner_interval),
          mean_min_f = mean(min_f_stat, na.rm = TRUE),
          method = "IVW",
          data_type = "Inner",
          .groups = "drop"
        ) %>%
        arrange(mean_mse)
    }
  }
}

#' Plot Simulation Results
#'
#' @param results Data frame from run_full_simulation_study()
#'
#' @return List of ggplot objects
plot_simulation_results <- function(results) {
  
  plots <- list()
  
  # 1. MSE comparison
  plots$mse <- results %>%
    filter(converged == TRUE) %>%
    ggplot(aes(x = scenario, y = poly_mse)) +
    geom_boxplot() +
    coord_flip() +
    labs(title = "Mean Squared Error by Scenario (Polynomial MVMR, Inner Interval)",
         x = "Scenario", y = "MSE") +
    theme_bw()
  
  # 2. Bias by quintile
  plots$bias_quintile <- results %>%
    filter(converged == TRUE) %>%
    select(scenario, starts_with("bias_q")) %>%
    pivot_longer(cols = starts_with("bias_q"),
                 names_to = "quintile",
                 values_to = "bias") %>%
    ggplot(aes(x = quintile, y = bias, color = scenario)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Bias Across Exposure Quintiles",
         x = "Quintile", y = "Bias") +
    theme_bw()
  
  # 3. Degree selection
  plots$degree_selection <- results %>%
    filter(converged == TRUE) %>%
    ggplot(aes(x = factor(selected_degree))) +
    geom_bar() +
    facet_wrap(~scenario) +
    labs(title = "Selected Polynomial Degree",
         x = "Degree", y = "Count") +
    theme_bw()
  
  # 4. R² performance
  plots$r_squared <- results %>%
    filter(converged == TRUE) %>%
    ggplot(aes(x = scenario, y = poly_r_squared)) +
    geom_boxplot() +
    coord_flip() +
    labs(title = "R² by Scenario (Polynomial MVMR, Inner Interval)",
         x = "Scenario", y = "R²") +
    theme_bw()
  
  # 5. Inner vs Full comparison (Polynomial)
  plots$inner_vs_full <- results %>%
    filter(converged == TRUE) %>%
    select(scenario, replicate, poly_mse, poly_mse_full) %>%
    pivot_longer(cols = c(poly_mse, poly_mse_full),
                 names_to = "data_type",
                 values_to = "mse") %>%
    mutate(data_type = ifelse(data_type == "poly_mse", "Inner Interval", "Full Data")) %>%
    ggplot(aes(x = scenario, y = mse, fill = data_type)) +
    geom_boxplot() +
    coord_flip() +
    scale_fill_manual(values = c("Inner Interval" = "steelblue", "Full Data" = "coral")) +
    labs(title = "MSE: Inner Interval vs Full Data (Polynomial MVMR)",
         x = "Scenario", y = "MSE", fill = "Data Type") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  # 6. Polynomial vs IVW comparison (Inner Interval)
  plots$poly_vs_ivw <- results %>%
    filter(converged == TRUE) %>%
    select(scenario, replicate, poly_mse, ivw_mse) %>%
    pivot_longer(cols = c(poly_mse, ivw_mse),
                 names_to = "method",
                 values_to = "mse") %>%
    mutate(method = ifelse(method == "poly_mse", "Polynomial MVMR", "IVW")) %>%
    ggplot(aes(x = scenario, y = mse, fill = method)) +
    geom_boxplot() +
    coord_flip() +
    scale_fill_manual(values = c("Polynomial MVMR" = "darkgreen", "IVW" = "darkorange")) +
    labs(title = "MSE: Polynomial MVMR vs Standard IVW (Inner Interval)",
         x = "Scenario", y = "MSE", fill = "Method") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  # 7. Method comparison by scenario
  plots$method_comparison_facet <- results %>%
    filter(converged == TRUE) %>%
    select(scenario, replicate, poly_mse, ivw_mse) %>%
    pivot_longer(cols = c(poly_mse, ivw_mse),
                 names_to = "method",
                 values_to = "mse") %>%
    mutate(method = ifelse(method == "poly_mse", "Polynomial MVMR", "IVW")) %>%
    ggplot(aes(x = method, y = mse, fill = method)) +
    geom_boxplot() +
    facet_wrap(~scenario, scales = "free_y") +
    scale_fill_manual(values = c("Polynomial MVMR" = "darkgreen", "IVW" = "darkorange")) +
    labs(title = "MSE Comparison by Scenario",
         x = "Method", y = "MSE") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(plots)
}

# Check for required source file
if (is.null(opt$source_file)) {
  stop("Error: --source_file is required. Please provide the path to your simulation functions file.")
}

if (!file.exists(opt$source_file)) {
  stop(paste("Error: Source file not found:", opt$source_file))
}

# Source the functions file
cat("Sourcing functions from:", opt$source_file, "\n")
source(opt$source_file)

# Parse comma-separated values into numeric vectors
prop_int <- as.numeric(strsplit(opt$prop_int, ",")[[1]])
inst_int <- as.integer(strsplit(opt$inst_int, ",")[[1]])
ss_int <- as.integer(strsplit(opt$ss_int, ",")[[1]])

# Print configuration
cat("=== Simulation Configuration ===\n")
cat("Inner interval:", opt$inner_interval, "\n")
cat("N replicates:", opt$n_replicates, "\n")
cat("Pleiotropy proportions:", paste(prop_int, collapse=", "), "\n")
cat("Instruments:", paste(inst_int, collapse=", "), "\n")
cat("Sample sizes:", paste(ss_int, collapse=", "), "\n")
cat("Output directory:", opt$output_dir, "\n")
cat("Output filename:", opt$output_name, "\n")
cat("================================\n\n")

# Create output directory if it doesn't exist
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
  cat("Created output directory:", opt$output_dir, "\n")
}

# Run simulation
cat("Starting simulation...\n")
results <- run_full_simulation_study(
  inner_interval = opt$inner_interval,
  n_replicates = opt$n_replicates,
  prop_int = prop_int,
  inst_int = inst_int,
  ss_int = ss_int
)

# Construct full output path
output_path <- file.path(opt$output_dir, opt$output_name)

# Save results
saveRDS(results, file = output_path)
cat("\nSimulation complete!\n")
cat("Results saved to:", output_path, "\n")
cat("Total rows:", nrow(results), "\n")

# Print summary
if (nrow(results) > 0) {
  cat("\n=== Quick Summary ===\n")
  cat("Converged simulations:", sum(results$converged, na.rm=TRUE), "\n")
  cat("Failed simulations:", sum(!results$converged, na.rm=TRUE), "\n")
}