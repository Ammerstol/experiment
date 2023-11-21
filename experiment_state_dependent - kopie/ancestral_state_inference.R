library(diversitree)
library(ape)

get_type <- function(name) {
  name_vec <- strsplit(name, " ")[[1]]
  return(as.numeric(name_vec[length(name_vec)]))
}

determine_parent_type <- function(parent_index, probas) { #parent index from sorted unique parents list
  return(probas[parent_index] <= 0.5)
}

filter_names <- function(all_names, params) {
  filtered_names <- c()
  for (name in all_names) {
    if (substring(name, 1, 1) %in% params) {
      filtered_names <- append(filtered_names, name)
    }
  }
  return(filtered_names)
}


run_trait_evolution_mcmc <- function(chain_length, matrix_i, exp, method, prior_mean = 2, seed = 1) {
  phy <- read.nexus(paste("experiment_real/real_tree_", matrix_i, "_exp_", exp, ".nex", sep = ""))
  
  tips <- phy$tip.label
  tip_df <- phy$edge[phy$edge[,2]<=length(tips),]
  n_tips <- nrow(tip_df)
  tip_parents_unicity_sorted <- sort(tip_df[,1][!duplicated(tip_df[,1])])
  
  parent_type <- integer(n_tips)
  child_type <- integer(n_tips)
  marginals <- numeric(n_tips)


  types <- c()
  for (name in tips) {
    types <- append(types, get_type(name))
  }
  states <- setNames(types, tips)
  k = max(types) + 1
  set.seed(seed)
  
  if (method == "mkn") {

    lik <- make.mkn(phy, states + 1, k)
    pars <- .2*rep(1, k*(k-1))
  } else if (method == "MuSSE") {
    lik <- make.musse(phy, states + 1, k)
    # mu_names <- argnames(lik)[(k+1):(k+k)]
    # for (mu_name in mu_names) { # restrict all extinction rates to be the same
    #   lik <- constrain(lik, formula(paste(mu_name, "~ mu")), extra="mu") 
    # }
    # pars <- .2*rep(1, k^2 + 1)
    pars <- .2*rep(1, k^2 + k)
  }

  prior <- make.prior.exponential(prior_mean)
  samples <- mcmc(lik, pars, chain_length, w=1, prior=prior,
                    print.every=0)
  print(samples)                  
  st_mcmc <- apply(samples[2:(length(pars) + 1)], 1, function(x) asr.marginal(lik, x, nodes = tip_parents_unicity_sorted - n_tips)[1,])
  
  if (method == "mkn") {  
    header_names <- filter_names(names(samples), c("l", "q"))
    } 
  else if (method == "MuSSE") {
    header_names <- filter_names(names(samples), c("m", "l", "q"))
  }
  params_inf <- colMeans(samples[header_names])
  return(params_inf)
}

mcmc_length <- 1000
experiments_and_matrix_indices <- read.delim("experiments_and_matrix_indices.txt", header = FALSE, sep = " ")
experiments_and_matrix_indices <- experiments_and_matrix_indices[(nrow(experiments_and_matrix_indices)-1),]
methods = c("mkn", "MuSSE")
method = methods[2]

params_inf <- run_trait_evolution_mcmc(1, experiments_and_matrix_indices[1,2], experiments_and_matrix_indices[1,1], method, seed = 1)
columns <- names(params_inf)
df <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(df) = columns
max_exp <- max(experiments_and_matrix_indices[, 1])

start <- Sys.time()
for (i in 1:nrow(experiments_and_matrix_indices)) {
  print(paste(i, " / ", nrow(experiments_and_matrix_indices)))
  exp = experiments_and_matrix_indices[i, 1]
  matrix_i = experiments_and_matrix_indices[i, 2]
  params_inf <- run_trait_evolution_mcmc(mcmc_length, matrix_i, exp, method, seed = max_exp*matrix_i + exp)
  df <- rbind(df,as.list(params_inf))
}


if (method == "mkn") {
  df$exp <- as.integer(experiments_and_matrix_indices[,1])
  df$matrix_i <- as.integer(experiments_and_matrix_indices[,2])
  write.csv(df, "experiment_real/output_real_tree_mkn.csv", row.names=FALSE)
} else if (method == "MuSSE") {
  df$exp <- as.integer(experiments_and_matrix_indices[,1])
  df$matrix_i <- as.integer(experiments_and_matrix_indices[,2])
  write.csv(df, "experiment_real/output_real_tree_MuSSE.csv", row.names=FALSE)
}