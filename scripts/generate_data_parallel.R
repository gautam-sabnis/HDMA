libs <- c("dplyr", "argparse", "tidyr", "mvtnorm", "Matrix")
sapply(libs, require, character.only = TRUE)

# Make folders where output will be stored.
setwd("/projects/kumar-lab/sabnig/HDMA/mediation_DNAm/")
dir.create("simulation_datasets")
dir.create("simulation_results")

#load variance-covariance matrix, which is stored as a sparse matrix with the lower triangle set to zero
load("vcv_sparse.rda") #loads object "vcv_sparse"
vcv_sparse <- as.matrix(vcv_sparse)
vcv <- vcv_sparse
vcv[lower.tri(vcv)] <- t(vcv)[lower.tri(vcv)]
rm(vcv_sparse)

parser <- argparse::ArgumentParser()
parser$add_argument("--seed", "-rng", type = "numeric")
parser$add_argument("--p", "-P", type = "numeric", help = "Number of features/variables")
parser$add_argument("--n", "-N", type = "numeric", help = "Number of samples")
parser$add_argument("--ndat", "-rep", type = "numeric", help = "Number of datasets")
parser$add_argument("--array_id", "-id", type = "numeric", help = "Array ID")

args <- parser$parse_args()
n <- args$n
p <- args$p
ndat <- args$ndat
seed <- args$seed
seed1 <- 35 #chosen to make residual variance positive
array_id <- args$array_id

tab <- tidyr::expand_grid(setting = c("bl", "hc", "ns"), dataset = c(1), seed2 = 1:ndat)

#array_id <- Sys.getenv('SLURM_ARRAY_TASK_ID') |> as.numeric()
setting <- tab[array_id, "setting"] |> pull()
d <- tab[array_id, "dataset"] |> pull()
seed2 <- tab[array_id, "seed2"] |> pull()

if (setting %in% c("bl", "hc")){

  seed <- seed

  na <- 80 #number of active features for alpha
  nb <- 80 #number of active features for beta
  nab <- 20 #number of active features for alpha and beta

  pve_a1 <- 0.2 
  pve_de1 <- 0.1
  pve_ie1 <- 0.1

  va <- 1
  valpha <- 1
  vbeta <- 1
  vm <- 1 #variance of mediators with alpha = 0. Helps to make this small.
  r <- 1 #parameter to regularize variance-covariance
    
  if(setting == "hc"){
    r <- 0.1
  }

  pve_a <- pve_a1
  pve_de <- pve_de1
  pve_ie <- pve_ie1
    
  if(d == 2){
    pve_a <- pve_a / 2
  } else if(d == 3){
    pve_de <- pve_de / 2
  } else if(d == 4){
    pve_ie <- pve_ie / 2
  }

  #SAMPLE POPULATION ATTRIBUTES
  set.seed(seed)
    
  #Sample alpha
  alpha_a <- rep(0,p)
  which_a <- sort(sample(p,na))
  alpha_a[which_a] <- as.numeric(scale(rnorm(na))) * sqrt(valpha)
    
  #Sample beta
  beta_m <- rep(0,p)
  which_ab <- sort(sample(which_a,nab))
  which_b <- sort(c(which_ab, sample(setdiff(1:p,which_a),nb - nab)))
  beta_m[which_b] <- as.numeric(scale(rnorm(nb))) * sqrt(vbeta)
  #beta_m[which_b] <- scale(c(round(rnorm(10, 8, 3), 2)*sample(c(-1,1), replace = TRUE, 10), rnorm(nb - 10))) * sqrt(vbeta)
      
  ab <- alpha_a * beta_m
  tie <- sum(ab)
    
  #Load mediator variance-covariance
  # vcv <- bigreadr::fread2("vcov2000.csv")
  set.seed(seed)
  which_mediators <- (sample(2000,p))
  vcv1 <- as.matrix(vcv[which_mediators,which_mediators])
    
  #standardize mediators
  sds <- sqrt(diag(vcv1))
  vcv1 <- t(vcv1 / sds) / sds
    
  #Regularize variance-covariance matrix so that it's not singular
  diag(vcv1) <- diag(vcv1) + r
    
  #Set desired mediator variance
  scale <- sqrt(diag(vcv1)) / sqrt(vm)
  vcv1 <- t(vcv1 / scale) / scale
    
  #Scale certain mediators to get desired pve_a
  scale1 <- rep(1,p)
  scale1[which_a] <- sqrt( ((va * alpha_a[which_a] ^ 2) / pve_a -
                                va * alpha_a[which_a]^2)/ vm)
  vcv2 <- t(vcv1 * scale1) * scale1
    
  #Compute Y variance components
  var_ie <- va * sum(ab)^2
  var_y <- var_ie / pve_ie #desired variance of Y
  beta_a <- sqrt(var_y * pve_de / va) * sign(tie)
  vc2 <- 2 * beta_a * sum(ab) * va
  vc4 <- sum( (beta_m %*% t(beta_m)) * vcv2)
  var_res <- var_y - va * beta_a^2 - vc2 - var_ie - vc4; var_res
    
  if(var_res <= 0){
    message("RESIDUAL VARIANCE NEGATIVE. CANCELING.")
    break
  }
    
  #SAMPLE DATASET
  A <- as.numeric(scale(rnorm(n))) * sqrt(va)

}

if(setting == "ns"){
    
  #ASSIGN PARAMETERS
  seed <- seed
  seed1 <- seed1 #chosen to make residual variance positive
  #p <- 2000
    
  na <- 80 #number of active features for alpha
  nb <- 80 #number of active features for beta
  nab <- 20 #number of active features for alpha and beta
    
  pve_a1 <- 0.2
  pve_de1 <- 0.1
  pve_ie1 <- 0.1
    
  va <- 1
  valpha <- 1
  vbeta <- 1
  vm <- 1 #variance of mediators with alpha = 0. Helps to make this small.
  r <- 1 #parameter to regularize variance-covariance
    
  vbeta0 <- 0.2^2 #difference needs to be small enough for BSLMM to struggle
  valpha0 <- 0.2^2
    
  pve_a <- pve_a1
  pve_de <- pve_de1
  pve_ie <- pve_ie1
    
  if(d == 2){
    pve_a <- pve_a / 2
      
  } else if(d == 3){
    pve_de <- pve_de / 2
      
  } else if(d == 4){
    pve_ie <- pve_ie / 2
      
  }
    
  #SAMPLE POPULATION ATTRIBUTES
    
  set.seed(seed)
  #Sample alpha
  alpha_a <- rep(0,p)
  which_a <- sort(sample(p,na))
  alpha_a[which_a] <- as.numeric(scale(rnorm(na))) * sqrt(valpha)
    
  #Sample beta
  beta_m <- rep(0,p)
  which_ab <- sort(sample(which_a,nab))
  which_b <- sort(c(which_ab, sample(setdiff(1:p,which_a),nb - nab)))
  beta_m[which_b] <- as.numeric(scale(rnorm(nb))) * sqrt(vbeta)
  #beta_m[which_b] <- scale(c(round(rnorm(10, 8, 3), 2)*sample(c(-1,1), replace = TRUE, 10), rnorm(nb - 10))) * sqrt(vbeta)
      
  #Fill in smaller effects
  set.seed(seed1)
  alpha_a[-which_a] <- as.numeric(scale(rnorm(p - na))) * sqrt(valpha0)
  beta_m[-which_b] <- as.numeric(scale(rnorm(p - nb))) * sqrt(vbeta0)
    
  ab <- alpha_a * beta_m
  tie <- sum(ab); tie
    
  # }
    
  sum(alpha_a[which_ab] * beta_m[which_ab])
    
  #Load mediator variance-covariance
  # vcv <- bigreadr::fread2("vcov2000.csv")
  set.seed(seed)
  which_mediators <- (sample(2000,p))
  vcv1 <- as.matrix(vcv[which_mediators,which_mediators])
    
  #standardize mediators
  sds <- sqrt(diag(vcv1))
  vcv1 <- t(vcv1 / sds) / sds
    
  #Regularize variance-covariance matrix so that it's not singular
  diag(vcv1) <- diag(vcv1) + r
    
  #Set desired mediator variance
  scale <- sqrt(diag(vcv1)) / sqrt(vm)
  vcv1 <- t(vcv1 / scale) / scale
    
  #Scale certain mediators to get desired pve_a
  scale1 <- sqrt( ((va * alpha_a ^ 2) / pve_a - va * alpha_a^2)/ vm)
  vcv2 <- t(vcv1 * scale1) * scale1
    
  #Compute Y variance components
  var_ie <- va * sum(ab)^2
  var_y <- var_ie / pve_ie #desired variance of Y
  beta_a <- sqrt(var_y * pve_de / va) * sign(tie)
  vc2 <- 2 * beta_a * sum(ab) * va
  vc4 <- sum( (beta_m %*% t(beta_m)) * vcv2)
  var_res <- var_y - va * beta_a^2 - vc2 - var_ie - vc4; var_res
    
  if(var_res <= 0){
    message("RESIDUAL VARIANCE NEGATIVE. CANCELING.")
    break
  }
    
  #SAMPLE DATASET
  A <- as.numeric(scale(rnorm(n))) * sqrt(va)
}

set.seed(seed2)
EM <- rmvnorm(n,rep(0,p),vcv2)
M <- A %*% t(alpha_a) + EM
    
ey <- rnorm(n,0,sqrt(var_res))
Y <- A * beta_a + M %*% beta_m + ey

save(Y, M, A, alpha_a, beta_m, beta_a, pve_a, pve_de, pve_ie,which_a, which_ab,which_b, file = paste0("simulation_datasets/sim_data_", setting,"_d",d,"_s",seed2,".rda"))
