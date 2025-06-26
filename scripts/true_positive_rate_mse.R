# Script to process global mediation effects
library(tidyverse)
library(bama)

rm(list = ls())
setwd("/projects/kumar-lab/sabnig/HDMA/mediation_DNAm/")

ndat <- 100 # Again, set this 100 if replicating the full study. This is equal to number of replicates parameter in generate_data_parallel.sh

# All settings to loop through
#settings <- c("bl", "hc", "ns","ln","ln_hc","ln_ns")
settings <- c("bl", "hc", "ns")
datasets <- c(1)
seed2 = 1:ndat

# Empty lists of results
mylist_tpr <- list()
mylist_mse0 <- list()
mylist_mse1 <- list()

# Function determine p-value cutoff for empirical FDR correction
test_fdr <- 
  function(score, truth = as.numeric(alpha_a * beta_m != 0), fdr = 0.1, cutoff.lower = T) {
    cutoffs <- sort(unique(score))
    
    
    
    if (cutoff.lower) {
      for (c in rev(cutoffs)) {
        sig <- score <= c
        
        fdp <- sum(sig & !truth) / sum(sig)
        
        if (fdp <= fdr) {
          break
        }
        
      }
      
    } else{
      for (c in cutoffs) {
        sig <- score >= c
        
        fdp <- sum(sig & !truth) / sum(sig)
        
        if (fdp <= fdr) {
          break
        }
        
      }
    }
    
    
    if(fdp <= fdr){
      return(sig)
      
    }else{ #this is for when the loop expires without finding a solution
      return(rep(F,length(sig)))
    }
    
    
    
  }

# Function to determine false discovery proportion for a given
# set of pathway LASSO estimates, since we have estimates for 
# every lambda. We choose the lowest lambda (i.e., the least
# shrunken model) attaining an FDP < 0.10
plasso_fdp <- function(estimates, which_cols, truth){
  # estimates: contribution estimates for a subset of mediators, some zero
  # which_cols: which_entries in "truth" were estimated at all 
  # truth: whether each of the 2000 mediators is active
  selected <- truth & F
  selected[which_cols] <- estimates != 0
  if(!any(selected)){return(1)} # Doesn't really matter what you return here.
  fdp <- sum(selected & !truth) / sum(selected)
  return(fdp)
  
}

setting <- 'hc'

for (setting in settings) {
  mylist_tpr[[setting]] <- list()
  mylist_mse0[[setting]] <- list()
  mylist_mse1[[setting]] <- list()
  
  for (dataset in datasets) {
    mylist_tpr[[setting]][[dataset]] <- list()
    mylist_mse0[[setting]][[dataset]] <- list()
    mylist_mse1[[setting]][[dataset]] <- list()
    
    for (seed in seed2) {
      temp <- ls()
      
      results <- list()
      
      # Load dataset
      load(paste0("simulation_datasets/sim_data_",setting,"_d",dataset,"_s",seed,".rda"))
      
      # Make empty dataset of results
      dat <- tibble(mediator = paste0("M", 1:2000),
                    beta_m,
                    alpha_a,
                    ab = alpha_a * beta_m)
      
      # Determine which mediators are active
      if (grepl("ns",setting)) {
        truth <- as.numeric((1:ncol(M)) %in% which_ab)
        
      } else{
        truth <- as.numeric(alpha_a * beta_m != 0)
      }
      
      # Load one-at-a-time results
      load(
        file = paste0(
          "simulation_results/one-at-a-time/sim_out_uni_",
          setting,
          "_d",
          dataset,
          "_s",
          seed,
          ".rda"
        )
      )
      
      results$uni <-
        out_uni |>
        rename(pv = ab_pv) |> 
        select(mediator,beta_hat,alpha_hat,ab_hat, pv) |> 
        left_join(x = dat, by = "mediator") |>
        mutate(method = "Univariate")
            
      # Load BSLMM results
      bslmm <- read_csv(paste0("simulation_results/bslmm/sim_out_bslmm_",setting,"_d", dataset,"_s",seed,".csv"))[,-1]
      results$bslmm <-
        bslmm |> 
        select(beta_hat, alpha_hat, pip) |> 
        mutate(
          ab_hat = alpha_hat * beta_hat, mediator = paste0("M",1:2000)) |>
          left_join(x = dat, by = "mediator") |>
          mutate(method = "BSLMM")

      #Load BOOM results
      boom <- read_csv(paste0("simulation_results/boom/sim_out_boom_",setting,"_d", dataset,"_s",seed,".csv"))[,-1]

      results$boom <-
        boom |> 
        select(beta_hat, alpha_hat, pip) |> 
        mutate(
          ab_hat = alpha_hat * beta_hat, mediator = paste0("M",1:2000)) |>
          left_join(x = dat, by = "mediator") |>
          mutate(method = "BOOM")
      
      d.mse0 <-
        results |>
        bind_rows() |>
        mutate(error = (ab_hat - ab) ^ 2, #error for MSE
               prb = 100 * abs(ab_hat - ab) / ab, #percent relative bias
               method = factor(
                 method,
                 levels = c('Univariate', 'BSLMM', 'BOOM')
)) |>
        group_by(method) |>
        filter(!truth) |>
        summarize(
          mse = mean(error)
        ) |>
        ungroup()  |> 
        mutate(
          rmse = mse / mse[method == "Univariate"],
          setting = setting,
          dataset = paste0("Setting ",dataset),
          seed2 = seed
        )
      
      d.mse1 <-
        results |>
        bind_rows() |>
        mutate(error = (ab_hat - ab) ^ 2, #error for MSE
               prb = 100 * abs(ab_hat - ab) / ab, #percent relative bias
               method = factor(
                 method,
                 levels = c('Univariate', 'BSLMM', 'BOOM')
               )) |>
        group_by(method) |>
        filter(as.logical(truth)) |>
        summarize(
          mse = mean(error),
        ) |>
        ungroup()  |> 
        mutate(
          rmse = mse / mse[method == "Univariate"],
          setting = setting,
          dataset = paste0("Setting ",dataset),
          seed2 = seed
        )
      
      d.tpr <-
        results |>
        bind_rows() |>
        mutate(score = ifelse(is.na(pv), 1 - pip, pv),
               method = factor(
                 method,
                 levels = c('Univariate', 'BSLMM', 'BOOM')
               )) |>
        group_by(method) |>
        mutate(sig = test_fdr(score, truth)) |>
        summarize(
          np = sum(sig),
          ntp = sum(sig & truth),
          nt = sum(truth),
          tpr = sum(sig & truth) / sum(truth),
          nfp = sum(sig & !truth),
          fpr = nfp / sum(!truth)
        ) |>
        ungroup() |> 
        mutate(setting = setting,
               dataset = paste0("Setting ",dataset),
               seed2 = seed)
      
      
      mylist_tpr[[setting]][[dataset]][[seed]] <- d.tpr
      mylist_mse0[[setting]][[dataset]][[seed]] <- d.mse0
      mylist_mse1[[setting]][[dataset]][[seed]] <- d.mse1
      
      #remove temporary objects (but is everything not temporary?)
      #rm(list = setdiff(ls(), temp))
      
    }
  }
}

# True positive rate
results |>
  rename(tpr1 = tpr) |> 
  group_by(setting, dataset, method) |>
  summarize(
    tpr = mean(tpr1),
    tpr_l = quantile(tpr1, 0.025),
    tpr_u = quantile(tpr1, 0.975)
  ) |>
  as.data.frame() |> 
  write_csv("simulation_results/tpr.csv")


# False positive rate - looked at this but not in the paper. They were low.
#results |>
#  rename(fpr1 = fpr) |> 
#  group_by(setting, dataset, method) |>
#  summarize(
#    fpr = mean(fpr1),
#    fpr_l = quantile(fpr1, 0.025),
#    fpr_u = quantile(fpr1, 0.975)
#  ) |>
#  as.data.frame() |> 
#  write_csv("simulation_results/fpr.csv")


# MSE for inactive mediators
#mylist_mse0 |>
#  map(bind_rows) |>
#  bind_rows() |>
#  rename(mse0 = mse, rmse0 = rmse) |> 
#  group_by(setting,dataset,method) |> 
#  summarize(
#    mse = mean(mse0),
#    mse_l = quantile(mse0, 0.025),
#    mse_u = quantile(mse0, 0.975),
#    rmse = mean(rmse0),
#    rmse_l = quantile(rmse0, 0.025),
#    rmse_u = quantile(rmse0, 0.975)
#  ) |>
#  as.data.frame() |> 
#  write_csv("simulation_results/mse_inactive.csv")


# MSE among active mediators
#mylist_mse1 |>
#  map(bind_rows) |>
#  bind_rows() |>
#  rename(mse0 = mse,rmse0 = rmse) |> 
#  group_by(setting,dataset,method) |> 
#  summarize(
#    mse = mean(mse0),
#    mse_l = quantile(mse0, 0.025),
#    mse_u = quantile(mse0, 0.975),
#    rmse = mean(rmse0),
#    rmse_l = quantile(rmse0, 0.025),
#    rmse_u = quantile(rmse0, 0.975)
#  ) |>
#  as.data.frame() |> 
#  write_csv("simulation_results/mse_active.csv")

tmp <- do.call(rbind, mylist_tpr$ns[[1]])

tmp |> ggplot(aes(x = method, y = tpr)) + geom_bar(stat = 'identity') + theme_classic() + labs(y = "True Positive Rate", x = "Method") 

tmp |> ggplot(aes(x = method, y = tpr)) + stat_summary(fun = mean, geom = 'bar') + stat_summary(fun.data = mean_se, geom = 'errorbar', width = 0.2)+ labs(y = "True Positive Rate", x = "") + theme_minimal(base_size = 14) + scale_fill_brewer(palette = "Set1") + expand_limits(x = 0, y = 0)
ggsave('simulation_results/plots/tpr_ns.pdf')

tmp <- do.call(rbind, mylist_mse0$ns[[1]])

tmp |> ggplot(aes(x = method, y = mse)) + geom_bar(stat = 'identity') + theme_classic() + labs(y = "MSE", x = "Method")

tmp |> ggplot(aes(x = method, y = mse)) + stat_summary(fun = mean, geom = 'bar') + stat_summary(fun.data = mean_se, geom = 'errorbar', width = 0.2)+ labs(y = "MSE", x = "") + theme_minimal(base_size = 14) + scale_fill_brewer(palette = "Set1") + expand_limits(x = 0, y = 0)
ggsave('simulation_results/plots/mse_ns.pdf')


tmp <- do.call(rbind, mylist_mse1$hc[[1]])

tmp |> ggplot(aes(x = method, y = mse)) + geom_bar(stat = 'identity') + theme_classic() + labs(y = "MSE", x = "Method")

tmp |> ggplot(aes(x = method, y = mse)) + stat_summary(fun = mean, geom = 'bar') + stat_summary(fun.data = mean_se, geom = 'errorbar', width = 0.2)+ labs(y = "MSE", x = "") + theme_minimal(base_size = 14) + scale_fill_brewer(palette = "Set1") + expand_limits(x = 0, y = 0)
ggsave('simulation_results/plots/mse_true_coef_hc.pdf')