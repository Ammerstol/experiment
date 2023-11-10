library(diversitree)
library(ape)

getTypeNum <- function(name){
  n <- nchar(name)
  type <- substr(name, n, n)
  if (type == "A") {
    return(0)
  }
  else {
    return(1)
  }
}

determine_parent_type <- function(parent_index, probas){ #parent index from sorted unique parents list
  return(probas[parent_index] <= 0.5)
}

run_mk2_mcmc <- function(proba_x_100, chain_length, exp, prior_mean = 2, seed = 1){
  phy <- read.nexus(paste("experiment_real/real_tree_", proba_x_100, "_exp_", exp, ".nex", sep = ""))
  
  tips <- phy$tip.label
  tip_df <- phy$edge[phy$edge[,2]<=length(tips),]
  n_tips <- nrow(tip_df)
  tip_parents_unicity_sorted <- sort(tip_df[,1][!duplicated(tip_df[,1])])
  
  types <- c()
  for (name in tips) {
    types <- append(types,getTypeNum(name))
  }
  states <- setNames(types, tips)
  pars <- c(.1, .2)
  set.seed(seed)
  
  lik.m2 <- make.mk2(phy, states)
  lik.m <- constrain(lik.m2, q10 ~ q01)
  fit.m <- find.mle(lik.m, pars[1], method="subplex")
  st_mle.m <- asr.marginal(lik.m, coef(fit.m), nodes = tip_parents_unicity_sorted - n_tips)[1,]
  
  prior <- make.prior.exponential(prior_mean)
  samples.m <- mcmc(lik.m, coef(fit.m)[1], chain_length, w=1, prior=prior,
                    print.every=0)
  st_mcmc.m <- apply(samples.m[2], 1, function(x) asr.marginal(lik.m, x, nodes = tip_parents_unicity_sorted - n_tips)[1,])
  st_mcmc.m.avg <- rowMeans(st_mcmc.m)
  
  parent_type <- integer(n_tips)
  child_type <- integer(n_tips)
  marginal_proba <- numeric(n_tips)
  
  for (row in 1:n_tips){
    parent_index <- match(tip_df[row,1], tip_parents_unicity_sorted)
    parent_type[row] <- determine_parent_type(parent_index, st_mcmc.m.avg)
    child_type[row] <- states[tip_df[row,2]]
    marginal_proba[row] <- st_mcmc.m.avg[parent_index]
  }
  
  jump_vector <- abs(parent_type-child_type)
  proba_inf <- sum(jump_vector)/n_tips
  rate_inf <- colMeans(samples.m[2])
  sum_marginals = sum(marginal_proba)
  weighted_proba_01 <- sum(marginal_proba * child_type)/sum_marginals
  weighted_proba_10 <- sum((1-marginal_proba) * (1-child_type))/(n_tips - sum_marginals)
  
  params_inf <- list("proba" = proba_inf, "rate" = rate_inf, "weighted_proba_01" = weighted_proba_01, "weighted_proba_10" = weighted_proba_10)
  
  return(params_inf)
}

mcmc_length <- 10000
probas <- read.delim("true_probabilities.txt", header = FALSE)[,1]
experiments <- read.delim("experiments.txt", header = FALSE)[,1]
probas_inferred <- c()
rates_inferred <- c()
weighted_proba_01 <- c()
weighted_proba_10 <- c()
proba_df <- c()
exp_df <- c()


start <- Sys.time()
for (exp in experiments){
  print(exp)
  for (proba in probas){
    proba_x_100 = as.integer(proba*100)
    params_inf <- run_mk2_mcmc(proba_x_100, mcmc_length, exp, seed = as.integer(length(experiments)*proba_x_100 + exp))
    probas_inferred <- append(probas_inferred, params_inf$proba)
    rates_inferred <- append(rates_inferred, params_inf$rate)
    weighted_proba_01 <- append(weighted_proba_01, params_inf$weighted_proba_01)
    weighted_proba_10 <- append(weighted_proba_10, params_inf$weighted_proba_10)
    proba_df <- append(proba_df, proba)
    exp_df <- append(exp_df, exp)
  }
}
print( Sys.time() - start )

df <- data.frame(exp_df, proba_df, probas_inferred, rates_inferred, weighted_proba_01, weighted_proba_10)
write.csv(df, "experiment_real/output_real_tree_accurracy.csv", row.names=FALSE)
