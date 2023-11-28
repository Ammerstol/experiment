setwd(getSrcDirectory(function(){})[1])

library(diversitree)
library(ape)
library(stringi)

get_type_synthetic_data <- function(name, method) {
  name_vec <- strsplit(name, " ")[[1]]
  if (method == "MuSSE") {
    return(as.numeric(name_vec[length(name_vec)]) + 1)
  }
  else {
    return(as.numeric(name_vec[length(name_vec)]))
  }
}

get_type_real_data <- function(age, age_groups) {
  type <- nrow(age_groups)
  for (age_i in 1:nrow(age_groups)) {
    if (age >= age_groups[age_i, "low"] & age <= age_groups[age_i, "up"]) {
      type <- age_i
      break
    }
  }
  return(type)
}

determine_parent_type <- function(parent_index, probas) { # parent index from sorted unique parents list
  return(probas[parent_index] <= 0.5)
}

filter_column_names <- function(all_names, params) { # to filter which parameters to keep as output
  filtered_names <- c()
  for (name in all_names) {
    if (substring(name, 1, 1) %in% params) {
      filtered_names <- append(filtered_names, name)
    }
  }
  return(filtered_names)
}

run_trait_evolution_mcmc <- function(chain_length, matrix_i, exp, method, age_groups, prior_mean = 2, seed = 1, take_mean = TRUE) {
  phy <- read.nexus(paste("experiment_real/", "real_tree_", matrix_i, "_exp_", exp, ".nex", sep = ""))

  tips <- phy$tip.label

  # tip_df <- phy$edge[phy$edge[,2]<=length(tips),]
  # n_tips <- nrow(tip_df)
  # tip_parents_unicity_sorted <- sort(tip_df[,1][!duplicated(tip_df[,1])])
  # parent_type <- integer(n_tips)
  # child_type <- integer(n_tips)
  # marginals <- numeric(n_tips)

  if (method != "bd") {
    types <- c()
    for (name in tips) {
      types <- append(types, get_type_synthetic_data(name, method))
    }
    states <- setNames(types, tips)
    k = max(types)
    set.seed(seed)
  }

  if (method == "mkn") {
    lik <- make.mkn(phy, states, k)
    pars <- .2*rep(1, k*(k-1))
  } else if (method == "MuSSE") {
    lik <- make.musse(phy, states, k)
    mu_names <- argnames(lik)[(k+1):(k+k)]
    for (mu_name in mu_names) { # restrict all extinction rates to be the same
      lik <- constrain(lik, formula(paste(mu_name, "~ mu")), extra="mu") 
    }
    pars <- .2*rep(1, k^2 + 1)
    # pars <- .2*rep(1, k^2 + k)
  } else if (method == "bd") {
    lik <- make.bd(phy)
    pars = .2*rep(1, 2)
  }

  prior <- make.prior.exponential(prior_mean)
  samples <- mcmc(lik, pars, chain_length, w=1, prior=prior,
                    print.every=0)  
  # st_mcmc <- apply(samples[2:(length(pars) + 1)], 1, function(x) asr.marginal(lik, x, nodes = tip_parents_unicity_sorted - n_tips)[1,])
  
  if (method == "mkn") {  
    header_names <- filter_column_names(names(samples), c("l", "q"))
  } 
  else if (method == "MuSSE") {
    header_names <- filter_column_names(names(samples), c("m", "l", "q"))
  } 
  else if (method == "bd") {
    header_names <- filter_column_names(names(samples), c("m", "l"))
  }

  if (take_mean) {
    params_inf <- colMeans(samples[header_names])
  } else {
     params_inf <- samples[header_names]
  }
  return(params_inf)
}

mcmc_length = 10000

methods = c("mkn", "MuSSE", "bd")
method = methods[2]
take_mean = FALSE

experiments_and_matrix_indices <- read.delim("experiments_and_matrix_indices.txt", header = FALSE, sep = " ")[1, ]
# experiments_and_matrix_indices <- experiments_and_matrix_indices[1,]

params_inf <- run_trait_evolution_mcmc(1, experiments_and_matrix_indices[1,1], experiments_and_matrix_indices[1,2], method, seed = 1)
columns <- names(params_inf)
df <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(df) = columns
max_exp <- max(experiments_and_matrix_indices[, 1])

start <- Sys.time()


if (take_mean) {
  for (i in 1:nrow(experiments_and_matrix_indices)) {
    print(paste(i, " / ", nrow(experiments_and_matrix_indices)))
    matrix_i = experiments_and_matrix_indices[i, 1]
    exp = experiments_and_matrix_indices[i, 2]
    params_inf <- run_trait_evolution_mcmc(mcmc_length, matrix_i, exp, method, seed = max_exp*matrix_i + exp)
    df <- rbind(df,as.list(params_inf))
  }
} else {
  matrix_i = experiments_and_matrix_indices[1, 1]
  exp = experiments_and_matrix_indices[1, 2]
  params_inf <- run_trait_evolution_mcmc(mcmc_length, matrix_i, exp, method, seed = max_exp*matrix_i + exp, take_mean = take_mean)
  df <- rbind(df,as.list(params_inf))
}

print(Sys.time() - start)

if (take_mean) {
  write.csv(df, paste("experiment_real/output_real_tree_", method,".csv", row.names=FALSE))
} else {
   write.csv(df, paste("experiment_real/output_real_tree_", method,"_single.csv", sep = ""), row.names=FALSE)
}