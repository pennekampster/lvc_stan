# Notes -------------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(here)

library(reshape2)

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear Workspace

source("../lvc_stan/functions2.R")

seed <- 154815468 # seed also used for stan # 462528635
set.seed(seed)

run <- FALSE

stan_code <- "/Users/pennekampster/Documents/Git projects/lvc_stan/Model/AdditiveErrorModel.stan" # Path to stan code
res_file <- "/Users/pennekampster/Documents/Git projects/lvc_stan/res_Chenyu.rds" # Path to results file (load or save)

n_it <- 2000
n_chains <- 3

# Data --------------------------------------------------------------------

# monoculture data
mono <- read.csv("../lvc_stan/Data/growth_curves_mono.csv")
mono <- mono %>% select(date, treatment, replicate, mean.dens.ml) %>% 
  mutate(NumDays = as.numeric(ymd(date)-ymd(unique(date)[1])+1),
         predict_spec = ifelse(treatment == "C", "Colp", "Dexio")) %>% rename(Age = NumDays, Volume = mean.dens.ml) %>% select(-date)

# polyculture data
load("../lvc_stan/Data/replicate_mean.Rdata")
poly <- replicate_mean
poly <- poly %>% rename(Age = NumDays, Volume = mean.dens.ml) %>% select(Age, treatment, replicate, Volume, predict_spec)

poly <- poly %>% dplyr::filter(!(Age==47 & treatment == "C+D" & predict_spec == "Dexio" & replicate == 1))

df <- rbind(mono, poly)

df<- df %>% filter(replicate %in% c(1,3) & predict_spec %in% c("Colp", "Dexio") & treatment %in% c("C+D", "C", "D")) %>% 
  mutate(Condition = ifelse(treatment == "D", "Single", ifelse(treatment == "C", "Single", "Mixed")),
         Experiment = ifelse(replicate==1, "A", "B")) %>% rename(Species = predict_spec) %>% select(-treatment)



# dfA <- read.csv("Data/Gause-yeist-competition-exp-1.csv", comment.char="#")
# dfB <- read.csv("Data/Gause-yeist-competition-exp-2.csv", comment.char="#")
# df2 <- process_data(dfA, dfB)

plot_data(df)
# ggsave("data.jpg", width = 28, height = 21, units = "cm", dpi = 300)

# Model -------------------------------------------------------------------

format_data <- function(df){
  df1_A <- subset(df, Condition == "C" & replicate == 1)
  df1_B <- subset(df, Condition == "C" & replicate == 3)
  
  df2_A <- subset(df, Condition == "D" & replicate == 1)
  df2_B <- subset(df, Condition == "D" & replicate == 3)
  
  df12_A <- subset(df, Condition == "C+D" & replicate == 1)
  df12_B <- subset(df, Condition == "C+D" & replicate == 3)
  
  list(
    D = 2,
    
    N1_A = nrow(df1_A),
    t1_A = df1_A$Age,
    y1_A = df1_A$Volume,
    
    N1_B = nrow(df1_B),
    t1_B = df1_B$Age,
    y1_B = df1_B$Volume,
    
    N2_A = nrow(df2_A),
    t2_A = df2_A$Age,
    y2_A = df2_A$Volume,
    
    N2_B = nrow(df2_B),
    t2_B = df2_B$Age,
    y2_B = df2_B$Volume,
    
    N12_A = nrow(df12_A)/2,
    t12_A = unique(df12_A$Age),
    y12_A = cbind(subset(df, Condition == "C+D" & Species == "Colp" & replicate == 1, Volume),
                  subset(df, Condition == "C+D" & Species == "Dexio"& replicate == 1, Volume)),
    
    N12_B = nrow(df12_B)/2,
    t12_B = unique(df12_B$Age),
    y12_B = cbind(subset(df, Condition == "C+D" & Species == "Colp" & replicate == 3, Volume),
                  subset(df, Condition == "C+D" & Species == "Dexio" & replicate == 3, Volume)),
    
    N_rep = 50,
    t_rep = 1:50
  )
}

data_stan <- format_data(df)

param <- c("r", "alpha", "sigma",
           "sigma_f0", "sigma_alpha",
           "k_uni", "k_multi",
           "f0",
           "f1_A", "f1_B", "f2_A", "f2_B", "f1_A", "f12_B",
           "y1_rep", "y2_rep", "y12_rep",
           "log_lik")

library(rstan)
# Parallel computing
rstan_options(auto_write = TRUE)
options(mc.cores = 3)

if (run){
  fit <- stan(file = stan_code, data = data_stan, iter = n_it, chains = n_chains, pars = param, seed = seed,
              control = list(adapt_delta = 0.99, max_treedepth = 12))
  
  saveRDS(fit, res_file)
}else{
  fit <- readRDS(res_file)
}

#saveRDS(fit, res_file)

# Explore solution -----------------------------------------------------------------

pdf(here("alphas.pdf"))
key_param <- c("alpha")
pairs(fit, pars = key_param)
plot(fit, pars = key_param, plotfun = "trace")
plot(fit, pars = c(key_param))
dev.off()

print(fit)

# out <- extract(fit, permuted = TRUE, inc_warmup = FALSE,
#         include = TRUE)
# 
# # extract parameter vector:
# alpha1 <- apply(out$alpha[,,1], 2, median)
# alpha2 <- apply(out$alpha[,,2], 2, median)
# r <- apply(out$r, 2, median)
# parms <- c(r, -alpha1[1], -alpha2[1], -alpha1[2], -alpha2[2])
# initialN <- c(188,331)
# 
# out2 <- deSolve::ode(y=initialN, times=seq(1, 60, length=100), func=lv_interaction, parms=parms)
# matplot(out2[,1], out2[,-1], type="l",
#         xlab="time", ylab="N", col=c("black",  "blue"), lty=c(1,3), lwd=2)
# legend("topright", c("Colp", "Dexio"), col=c("black",  "blue"), lwd=2, lty=c(1,3))
# points(subset(df, Condition == "Mixed" & Species == "Colp")$Age, subset(df, Condition == "Mixed" & Species == "Colp")$Volume, col="blue")
# points(subset(df, Condition == "Mixed" & Species == "Dexio")$Age, subset(df, Condition == "Mixed" & Species == "Dexio")$Volume, col="black")
# 


shinystan::launch_shinystan(fit)

# Posterior Predictive Checks -------------------------------------------------------------------

# rep <- process_replications_density(data_stan, fit, maxVolume = 15)
rep <- process_replications_spaghetti(data_stan, fit, draws = 100)

pdf(here("dynamics.pdf"))
plot_PPC(rep, df)
dev.off()
# Compare models ----------------------------------------------------------

library(loo)

log_lik1 <- extract_log_lik(fit, merge_chains = FALSE)
r_eff1 <- relative_eff(exp(log_lik1), chain_id = rep(1:dim(log_lik1)[2], each = dim(log_lik1)[1]))
(loo1 <- loo(log_lik1, r_eff = r_eff1))
# plot(loo1)

fit2 <- readRDS("Results/fit_LV12_exp2_v3.rds")
log_lik2 <- extract_log_lik(fit2, merge_chains = FALSE)
r_eff2 <- relative_eff(exp(log_lik2), chain_id = rep(1:dim(log_lik2)[2], each = dim(log_lik2)[1]))
(loo2 <- loo(log_lik2, r_eff = r_eff2))
comp <- compare(loo1, loo2)
print(comp) # negative elpd favors first model, to compared with SE
