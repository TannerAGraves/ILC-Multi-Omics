library(ggplot2)
setwd("C:/Users/JIC/Documents/UNIPD/Stat Learning/proj")
omics <- read.table("data.csv", header = TRUE, sep=",")
rs.cols <- grep(paste0("^", "rs"), names(omics))
cn.cols <- grep(paste0("^", "cn"), names(omics))
mu.cols <- grep(paste0("^", "mu"), names(omics))
pp.cols <-grep(paste0("^", "pp"), names(omics))

# Missing value imputation
# due to high dim, removing missing not an option; every row and col has missing values
rs <- omics[rs.cols]
pp <- omics[pp.cols]
rs.missing <- sum(apply(omics[rs.cols], 1, function(x) any(x == 0))) 
rs.medians <- apply(rs, 2, function(x) median(x[x != 0]))
pp.medians <- apply(pp, 2, function(x) median(x[x != 0]))
# Missing value rate
rs.missing.rate <- sum(rs == 0) / (length(rs.cols) * length(rs))
rs.missing.rate
pp.missing.rate <- sum(pp == 0) / (length(pp.cols) * length(pp))
pp.missing.rate

for(i in seq_along(rs.medians)) {
  rs[rs[, i] == 0, i] <- rs.medians[i]
}
for(i in seq_along(pp.medians)) {
  pp[pp[, i] == 0, i] <- pp.medians[i]
}
omics[c(rs.cols,pp.cols)] <- c(rs,pp)

# Feature distributions
# Gene Expression
long_rs <- tidyr::gather(omics[rs.cols], value = "value", rs.cols)
ggplot(long_rs, aes(x = value)) + 
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Smoothed Distribution of Gene Expression Levels", x = "Level", y = "Density")
# Protein Abundance
long_pp <- tidyr::gather(omics[pp.cols], value = "value", pp.cols)
ggplot(long_pp, aes(x = value)) + 
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Smoothed Distribution of Protein Abundance Levels", x = "Level", y = "Density")
# Copy Number
long_cn <- tidyr::gather(omics[cn.cols], value = "value", cn.cols)$value
barplot(table(long_cn)/sum(table(long_cn)))
# Mutations
long_mu <- tidyr::gather(omics[mu.cols], value = "value", mu.cols)$value
barplot(table(long_mu)/sum(table(long_mu)))
# Vital Status (Target)
barplot(table(omics$vital.status)/sum(table(omics$vital.status)))

# Normalize Continuous Features (Expression and Protein Abundance)
omics_norm <- omics
omics_norm[c(rs.cols,pp.cols)] <- scale(omics[c(rs.cols,pp.cols)])
colMeans(omics)[1:5]

# Ridge Regression
library(glmnet)
X <- omics_norm[,-ncol(omics_norm)]
y <- omics_norm$vital.status 
glm.out <- glm(vital.status~., data=omics_norm, family=binomial)

grid <- 10^seq(10, -2, length=100)
lasso.mod <- glmnet(X,y,alpha=1)
plot(lasso.mod, xvar="lambda", label=TRUE)
plot(lasso.mod, xvar="norm", label=TRUE)
plot(lasso.mod, xvar="dev", label=TRUE)
