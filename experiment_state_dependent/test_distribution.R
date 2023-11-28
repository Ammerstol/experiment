setwd(getSrcDirectory(function(){})[1])
df <- read.csv("experiment_real/output_real_tree_bd_single.csv")

lambda <- df[["lambda"]]
log_lambda <- log(lambda)

print(shapiro.test(lambda[seq(5001,10000)]))
print(shapiro.test(log_lambda[seq(5001,10000)]))