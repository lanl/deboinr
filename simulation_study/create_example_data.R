# Create data example for deboinr using Justin's breakthrough curves.

########################
## Code from Justin to get the PDF data:
####
## Load data
####
# Load feature .csv
setwd("~/Documents/UNEs/code_from_justin/current_lqd_biascorrection")
source("custom_functions_lqd.R")
wsp_feat_rs = read.csv("wsp_feat_rs.csv", header=FALSE)

# Load breakthrough curve data
graph_bt = list()
dfn_bt = list()

# Murph 6/29/2023: I'm going to follow along and create the LQD descrepencies 
# for the graph features.
setwd("~/Documents/UNEs/code_from_justin/current_lqd_biascorrection/tpl")
idx = 1
tpl_files = list.files()
for(run_dir in tpl_files){
  if(run_dir=="tpl_x218"){
    next
  }
  setwd(run_dir)
  
  graph_bt[[idx]] = read.csv("graph_times.dat", header=FALSE)[,1]
  dfn_bt[[idx]] = read.csv("dfn_times.dat", header=FALSE)[,1]
  
  setwd("~/Documents/UNEs/code_from_justin/current_lqd_biascorrection/tpl")
  idx = idx + 1
  print(idx)
}

setwd("~/DeBoinR/R")

# Create full data list from original
data_list = list()
for(idx in 1:length(graph_bt)){
  data_list[[idx]] = list(TG=graph_bt[[idx]], TD=dfn_bt[[idx]], 
                          feat=as.numeric(wsp_feat_rs[idx,]))
}

#####
## Split into training and test sets
#####
set.seed(14)

# Training and test indices
p_train = 0.8
train_idx = sort(sample(1:length(data_list), floor(p_train*length(data_list)), 
                        replace=FALSE))
test_idx = setdiff(1:length(data_list), train_idx)
n_train = length(train_idx)
n_test = length(test_idx)

# Training and test dataframes
train_list = data_list[train_idx]
test_list = data_list[test_idx]

#####
## Estimate supports from training
#####
# Settings and initialization of vectors
kappa_L = 2500
kappa_U = 1

range_lTG = matrix(NA, nrow=2, ncol=n_train)
range_lTD = matrix(NA, nrow=2, ncol=n_train)
range_lfe = matrix(NA, nrow=2, ncol=n_train)
sd_lTG = rep(NA, n_train)
sd_lTD = rep(NA, n_train)
sd_lfe = rep(NA, n_train)
n_lTG = rep(NA, n_train)
n_lTD = rep(NA, n_train)
n_lfe = rep(NA, n_train)

for(idx in 1:n_train){
  range_lTG[,idx] = quantile(log(train_list[[idx]]$TG), probs=c(0,1))
  range_lTD[,idx] = quantile(log(train_list[[idx]]$TD), probs=c(0,1))
  range_lfe[,idx] = quantile(log(train_list[[idx]]$feat), probs=c(0,1))
  sd_lTG[idx] = sd(log(train_list[[idx]]$TG))
  sd_lTD[idx] = sd(log(train_list[[idx]]$TD))
  sd_lfe[idx] = sd(log(train_list[[idx]]$feat))
  n_lTG[idx] = length(train_list[[idx]]$TG)
  n_lTD[idx] = length(train_list[[idx]]$TD)
  n_lfe[idx] = length(train_list[[idx]]$feat)
}

# Approach 1: common support for both graph and DFN log breakthroughs
min_L = min(c(min(range_lTG[1,]), min(range_lTD[1,])))
max_U = max(c(max(range_lTG[2,]), max(range_lTD[2,]))) 
pooled_sd = sum(c((n_lTG-1)*(sd_lTG^2), (n_lTD-1)*(sd_lTD^2)))
pooled_sd = sqrt(pooled_sd/(sum(n_lTG)+sum(n_lTD)-2*n_train))
omega_L = min_L - kappa_L*pooled_sd/sqrt(sum(n_lTG)+sum(n_lTD))
omega_U = max_U + kappa_U*pooled_sd/sqrt(sum(n_lTG)+sum(n_lTD))
#omega_U = 18

# Approach 2: separate supports for graph and DFN log breakthroughs
min_LG = min(range_lTG[1,])
max_UG = max(range_lTG[2,])
pooled_sdG = sqrt(sum((n_lTG-1)*(sd_lTG^2))/(sum(n_lTG)-n_train))
omega_LG = min_LG - kappa_L*pooled_sdG/sqrt(sum(n_lTG))
omega_UG = max_UG + kappa_U*pooled_sdG/sqrt(sum(n_lTG))

min_LD = min(range_lTD[1,])
max_UD = max(range_lTD[2,])
pooled_sdD = sqrt(sum((n_lTD-1)*(sd_lTD^2))/(sum(n_lTD)-n_train))
omega_LD = min_LD - kappa_L*pooled_sdD/sqrt(sum(n_lTD))
omega_UD = max_UD + kappa_U*pooled_sdD/sqrt(sum(n_lTD))

# Murph 06/29/2023: I don't need the LQD on graph features to be in this pool.
min_Lfe = min(range_lfe[1,])
max_Ufe = max(range_lfe[2,])
pooled_sdfe = sqrt(sum((n_lfe-1)*(sd_lfe^2))/(sum(n_lfe)-n_train))
omega_Lfe = min_Lfe - kappa_L*pooled_sdD/sqrt(sum(n_lfe))
omega_Ufe = max_Ufe + kappa_U*pooled_sdD/sqrt(sum(n_lfe))

## Conclusion: not much difference, just use approach 1 here

#####
## Compute PDFs, CDFs, and log quantile densities (LQDs)
#####
# Settings
n = 200
x_full_std = seq(from=0, to=1, length.out=n)
bw = 0.03 # bandwidth for Gaussian kernel density estimator (pre-specified, fixed)
alpha = 0.25 # mixture weight for numerically stable inverting of LQD

# Training
pdf_data = NULL
for(idx in 1:n_train){
  # Standardize based on support (omega_L, omega_U)
  train_list[[idx]]$lTG_std = (log(train_list[[idx]]$TG)-omega_L)/(omega_U-omega_L)
  train_list[[idx]]$lTD_std = (log(train_list[[idx]]$TD)-omega_L)/(omega_U-omega_L)
  train_list[[idx]]$lfe_std = (log(train_list[[idx]]$fe)-omega_Lfe)/(omega_Ufe-omega_Lfe)
  
  ##########
  ## PDFs ##
  # Original (standardized support)
  train_list[[idx]]$fG           = bkde(train_list[[idx]]$lTG_std, kernel="normal", bandwidth=bw,
                                        gridsize=n, range.x=c(0,1))$y
  neg_idx                        = which(train_list[[idx]]$fG<0)
  train_list[[idx]]$fG[neg_idx]  = 0
  train_list[[idx]]$fD           = bkde(train_list[[idx]]$lTD_std, kernel="normal", bandwidth=bw,
                                        gridsize=n, range.x=c(0,1))$y
  neg_idx                        = which(train_list[[idx]]$fD<0)
  train_list[[idx]]$fD[neg_idx]  = 0
  pdf_data                       = rbind(pdf_data, train_list[[idx]]$fD)
  
  
  train_list[[idx]]$ffe          = bkde(train_list[[idx]]$lfe_std, kernel="normal", bandwidth=bw,
                                        gridsize=n, range.x=c(0,1))$y
  neg_idx                        = which(train_list[[idx]]$ffe<0)
  train_list[[idx]]$ffe[neg_idx] = 0
  #train_list[[idx]]$fG = density(train_list[[idx]]$lTG_std, bw=bw, kernel="gaussian", 
  #                               from=0, to=1, n=n)$y
  #train_list[[idx]]$fD = density(train_list[[idx]]$lTD_std, bw=bw, kernel="gaussian", 
  #                               from=0, to=1, n=n)$y
  
  # Mixture with uniform
  train_list[[idx]]$fG_star  = alpha_mix(train_list[[idx]]$fG, alpha)
  train_list[[idx]]$fD_star  = alpha_mix(train_list[[idx]]$fD, alpha)
  train_list[[idx]]$ffe_star = alpha_mix(train_list[[idx]]$ffe, alpha)
  
  ##########
  ## CDFs ##
  # Original (standardized support)
  train_list[[idx]]$FG  = pdf_to_cdf(train_list[[idx]]$fG, x_full_std, norm=TRUE)
  train_list[[idx]]$FD  = pdf_to_cdf(train_list[[idx]]$fD, x_full_std, norm=TRUE)
  train_list[[idx]]$Ffe = pdf_to_cdf(train_list[[idx]]$ffe, x_full_std, norm=TRUE)
  
  # Mixture with uniform
  train_list[[idx]]$FG_star  = pdf_to_cdf(train_list[[idx]]$fG_star, x_full_std, 
                                          norm=TRUE)
  train_list[[idx]]$FD_star  = pdf_to_cdf(train_list[[idx]]$fD_star, x_full_std, 
                                          norm=TRUE)
  train_list[[idx]]$Ffe_star = pdf_to_cdf(train_list[[idx]]$ffe_star, x_full_std, 
                                          norm=TRUE)
  
  ########################
  ## Quantile functions ##
  ## Taken as inverse of CDF via linear interpolation
  # Original (standardized support)
  train_list[[idx]]$QG  = cdf_to_quant(train_list[[idx]]$FG, x_full_std)
  train_list[[idx]]$QD  = cdf_to_quant(train_list[[idx]]$FD, x_full_std)
  train_list[[idx]]$Qfe = cdf_to_quant(train_list[[idx]]$Ffe, x_full_std)
  
  # Mixture with uniform
  train_list[[idx]]$QG_star  = cdf_to_quant(train_list[[idx]]$FG_star, x_full_std)
  train_list[[idx]]$QD_star  = cdf_to_quant(train_list[[idx]]$FD_star, x_full_std)
  train_list[[idx]]$Qfe_star = cdf_to_quant(train_list[[idx]]$Ffe_star, x_full_std)
  
  ##########
  ## LQDs ##
  # Original (standardized support)
  train_list[[idx]]$psiG  = pdf_to_lqd(train_list[[idx]]$fG, train_list[[idx]]$QG, 
                                       x_full_std)
  train_list[[idx]]$psiD  = pdf_to_lqd(train_list[[idx]]$fD, train_list[[idx]]$QD, 
                                       x_full_std)
  train_list[[idx]]$psife = pdf_to_lqd(train_list[[idx]]$ffe, train_list[[idx]]$Qfe, 
                                       x_full_std)
  
  # Mixture with uniform
  train_list[[idx]]$psiG_star  = pdf_to_lqd(train_list[[idx]]$fG_star, 
                                            train_list[[idx]]$QG_star, x_full_std)
  train_list[[idx]]$psiD_star  = pdf_to_lqd(train_list[[idx]]$fD_star, 
                                            train_list[[idx]]$QD_star, x_full_std)
  train_list[[idx]]$psife_star = pdf_to_lqd(train_list[[idx]]$ffe_star, 
                                            train_list[[idx]]$Qfe_star, x_full_std)
  
  ###################
  ## Discrepancies ##
  # Mixture with uniform
  train_list[[idx]]$delta_star = train_list[[idx]]$psiD_star-train_list[[idx]]$psiG_star
}
pdf_data = rbind(unlist(x_full_std), pdf_data)
write.csv(pdf_data, "pdf_example_data.csv")


