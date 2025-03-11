#Script for implementing SpikeandSlab prior on simulated data
#Note that, before running this, one must generate data using "generate_data_parallel.R".

libs <- c("dplyr", "BoomSpikeSlab", "tidyr", "argparse")
lapply(libs, require, character.only = TRUE)

setwd("/projects/kumar-lab/sabnig/HDMA/mediation_DNAm/")

parser <- argparse::ArgumentParser()
parser$add_argument("--ndat", "-rep", type = "numeric", help = "Number of datasets")
parser$add_argument("--array_id", "-id", type = "numeric", help = "Array ID")

args <- parser$parse_args()
ndat <- args$ndat
array_id <- args$array_id

tab <- tidyr::expand_grid(setting = c("bl", "hc", "ns"), dataset = c(1), seed2 = 1:ndat)

setting <- tab[array_id, "setting"] |> pull()
d <- tab[array_id, "dataset"] |> pull()
seed2 <- tab[array_id, "seed2"] |> pull()

load(paste0("simulation_datasets/sim_data_",setting,"_d",d,"_s",seed2,".rda"))

C <- matrix(1,nrow(M),1) # covariate matrix

y_sd <- sd(Y)
m_sds <- apply(M, 2, sd)
a_sd <- sd(A)  
M <- as.matrix(scale(M))
Y <- as.numeric(scale(Y))
A <- as.numeric(scale(A))

alpha_hats <- apply(M, 2, function(x){summary(lm(x ~ A))$coefficients[2,1]})
alpha_hats <- alpha_hats*m_sds/a_sd
alpha_pvals <- apply(M, 2, function(x){summary(lm(x ~ A))$coefficients[2,4]})
alpha_pvals <- p.adjust(alpha_pvals, method = "fdr")

trueind_alpha <- which(alpha_a != 0)
thresh <- 0.05 
estind_alpha <- which(alpha_pvals < thresh)

TP_alpha <- length(intersect(estind_alpha, trueind_alpha))
FP_alpha <- length(setdiff(estind_alpha, trueind_alpha))
FN_alpha <- length(setdiff(trueind_alpha, estind_alpha))

mse_alpha <- mean((alpha_a[trueind_alpha] - alpha_hats[trueind_alpha])^2)


data <- data.frame(Y = Y, A = A, M = M)
outcome_model <- lm.spike(Y ~ A + M, data = data, niter = 40000)
dfbetam <- SummarizeSpikeSlabCoefficients(outcome_model$beta, burn = 10000) |> data.frame()
dfbetam <- dfbetam[!rownames(dfbetam) %in% c("A", "(Intercept)"),]
dfbetam$Mediator <- gsub("M", "", rownames(dfbetam)) |> as.numeric()
dfbetam <- dfbetam |> arrange(Mediator) |> as_tibble()
dfbetam <- dfbetam |> mutate(mean.inc = mean.inc * y_sd / m_sds)


#Inference
trueind_betam <- which(beta_m != 0) #true non-zero indices
thres <- 0.50 #inclusion prob threshold
estind_betam <- dfbetam[dfbetam$inc.prob > thres, ] |> pull(Mediator)


#inferred non-zero indices
TP <- length(intersect(estind_betam, trueind_betam))
FP <- length(setdiff(estind_betam, trueind_betam))
FN <- length(setdiff(trueind_betam, estind_betam))


mse <- mean((beta_m[trueind_betam] - dfbetam[trueind_betam, ] |> pull('mean.inc'))^2) #mean squared error for beta_m
dfcoef <- data.frame(True = beta_m[trueind_betam], Estimated = dfbetam[trueind_betam, ] |> pull('mean.inc'))

results <- list(TP = TP, FP = FP, FN = FN, mse = mse, dfcoef = dfcoef, TP_alpha = TP_alpha, FP_alpha = FP_alpha, FN_alpha = FN_alpha, mse_alpha = mse_alpha)

#beta_hat <- dfbetam[!rownames(dfbetam) %in% c("A", "(Intercept)"), "mean.inc"] * y_sd / m_sds
beta_hat <- dfbetam |> pull(mean.inc)
alpha_hat <- alpha_hats
#pip <- dfbetam[!rownames(dfbetam) %in% c("A", "(Intercept)"), "inc.prob"]
pip <- dfbetam |> pull(inc.prob) * as.numeric(alpha_pvals < 0.05)
tie_hat <- beta_hat * alpha_hat

out <- data.frame(beta_hat = beta_hat, alpha_hat = alpha_hat, pip = pip, tie_hat = tie_hat)

if(!dir.exists("/projects/kumar-lab/sabnig/HDMA/mediation_DNAm/simulation_results/boom")){
    dir.create("/projects/kumar-lab/sabnig/HDMA/mediation_DNAm/simulation_results/boom")
}

#saveRDS(results, file = paste0("/projects/kumar-lab/sabnig/HDMA/mediation_DNAm/simulation_results/boom/sim_out_boom_", setting, "_d", d, "_s", seed2, ".rds"))

out |> write.csv(file = paste0("/projects/kumar-lab/sabnig/HDMA/mediation_DNAm/simulation_results/boom/sim_out_boom_", setting, "_d", d, "_s", seed2, ".csv"))
