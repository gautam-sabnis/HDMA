#Script for implementing BSLMM on simulated data
#Note that, before running this, one must have already generated data using "generate_data_parallel.R".

libs <- c("tidyr", "bama", "argparse", "dplyr")
sapply(libs, require, character.only = TRUE)

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
  
# Bslmm tends to work better when scaling, re-scaling
y_sd <- sd(Y)
m_sds <- apply(M, 2, sd)
a_sd <- sd(A)
M <- as.matrix(scale(M))
Y <- as.numeric(scale(Y))
A <- as.numeric(scale(A))

# Get one-at-a-time coefficients so that we can input a good prior variance
alpha_hats <- apply(M, 2, function(x){summary(lm(x ~ A))$coefficients[2,1]})
alpha_var1 <- var(alpha_hats[abs(alpha_hats) > quantile(abs(alpha_hats),0.9)])
beta_hats <- apply(M, 2, function(x){summary(lm(Y ~ x + A))$coefficients[2,1]})
beta_var1 <- var(beta_hats[abs(beta_hats) > quantile(abs(beta_hats),0.9)])

out_bslmm <- bama(
    Y = as.vector(Y),
    A = A,
    C1 = C,
    C2 = C,
    M = M,
    burnin = 15000,
    ndraws = 20000,
    method = "BSLMM",
    control = list(k = 2, lma1 = alpha_var1, lm1 = beta_var1,
    l = 1, lambda0 = 0.04, lambda1 = 0.2, lambda2 = 0.2, 
    phi0 = 0.01, phi1 = 0.01, a0 = 0.01 * ncol(M),
    a1 = 0.05 * ncol(M), a2 = 0.05 * ncol(M),
    a3 = 0.89 * ncol(M)),seed = 123
)
  
out <- data.frame(beta_hat = colMeans(out_bslmm$beta.m) * y_sd / m_sds, alpha_hat = colMeans(out_bslmm$alpha.a) * m_sds / a_sd,tie_hat =c(with(out_bslmm, mean(rowSums(alpha.a * beta.m))) * y_sd / a_sd, rep(NA,1999)), pip = colMeans(with(out_bslmm,r1 * r3)))
  
write.csv(out, file = paste0("simulation_results/bslmm/sim_out_bslmm_",setting,"_d",d,"_s",seed2,".csv"))
