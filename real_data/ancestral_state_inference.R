# setwd("C:/Users/sjoerd/Downloads/experiment-state_dependent/experiment_state_dependent")

library(diversitree)
library(ape)
library(stringi)

get_type_synthetic_data <- function(name) {
  name_vec <- strsplit(name, " ")[[1]]
  return(as.numeric(name_vec[length(name_vec)]))
}

get_type_real_data <- function(age, age_groups) {
  type <- nrow(age_groups)
  print(age)
  for (age_i in 1:nrow(age_groups)) {
    print(age >= age_groups[age_i, "low"] & age <= age_groups[age_i, "up"])
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

run_trait_evolution_mcmc <- function(chain_length, country, week, method, age_groups, prior_mean = 2, seed = 1) {
  phy <- read.nexus(paste("tree_data/", country, "_wk_", week, "_beast_tree.nex", sep = ""))
  meta_df <- read.delim(paste("gisaid_data/", country, "_wk_", week, "_metadata.tsv", sep = ""), sep = "\t")

  tips <- phy$tip.label

  # tip_df <- phy$edge[phy$edge[,2]<=length(tips),]
  # n_tips <- nrow(tip_df)
  # tip_parents_unicity_sorted <- sort(tip_df[,1][!duplicated(tip_df[,1])])
  # parent_type <- integer(n_tips)
  # child_type <- integer(n_tips)
  # marginals <- numeric(n_tips)

  types <- c()
  for (name in tips) {
    print(name)
    age = meta_df[meta_df$strain == name, "age"]
    types <- append(types, get_type_real_data(age, age_groups))
  }
  print(age_groups)
  states <- setNames(types, tips)
  k = max(types) + 1
  set.seed(seed)
  
  if (method == "mkn") {
    lik <- make.mkn(phy, states + 1, k)
    pars <- .2*rep(1, k*(k-1))
  } else if (method == "MuSSE") {
    lik <- make.musse(phy, states + 1, k)
    mu_names <- argnames(lik)[(k+1):(k+k)]
    for (mu_name in mu_names) { # restrict all extinction rates to be the same
      lik <- constrain(lik, formula(paste(mu_name, "~ mu")), extra="mu") 
    }
    pars <- .2*rep(1, k^2 + 1)
    # pars <- .2*rep(1, k^2 + k)
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
  params_inf <- colMeans(samples[header_names])
  return(params_inf)
}

age_groups = data.frame(low = seq(0, 90, 10), up = seq(9, 99, 10))
methods = c("mkn", "MuSSE")
countries = c("slovenia", "belgium", "luxembourg")
weeks = c("1", "2", "1_2")

method = methods[2]
country = countries[1]
week = weeks[2]


start <- Sys.time()

params_inf <- run_trait_evolution_mcmc(mcmc_length, country, week, method, age_groups, seed = sum(utf8ToInt(country))*sum(utf8ToInt(week)))
df <- rbind(df,as.list(params_inf))

print(Sys.time() - start)

if (method == "mkn") {
  write.csv(df, "experiment_real/output_real_tree_mkn.csv", row.names=FALSE)
} else if (method == "MuSSE") {
  write.csv(df, "experiment_real/output_real_tree_MuSSE.csv", row.names=FALSE)
}
