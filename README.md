# coinflibs
**Co**nditional **Inf**erence after **Li**kelihood-**b**ased **S**election

The package contains functions to calculate limits and conduct inference in a selective manner for linear models after likelihood- or test-based model selection. 

Please see https://arxiv.org/abs/1706.09796 for more details.

## Example
 
```R
# install and load package
install_github("davidruegamer/coinflibs")
library("coinflibs")

library(MASS)
data("cpus")
# Fit initial model
cpus$perf <- log10(cpus$perf)
cpus$cach <- as.factor(cpus$cach)
mod <- lm(perf ~ .-name, data = cpus)

# use the stepAIC function to find the best model in a backward 
# stepwise search
cpus.lm <- stepAIC(mod, trace = FALSE, direction = "backward", steps = 3)
# check model selection
cpus.lm$anova$Step

# recalculate all visited models in the first step
lom1 <- c(lapply(attr(mod$terms, "term.labels"), function(x) 
update(mod, as.formula(paste0("perf ~ .-", x)))), list(mod))

# perform likelihood ratio test at level 
alpha = 0.001

# check for non-significant variables
coefTable <- anova(cpus.lm)
drop <- rownames(coefTable)[alpha < coefTable[-nrow(coefTable),5]]

# drop non-significant variable
cpus.lm2 <- update(cpus.lm, as.formula(paste0(".~.-",drop)))

# compute selective inference
selinf(list(cpus.lm, cpus.lm2), lom1, 
response = cpus$perf,
what = c("Ftest", "aic"),
sd = summary(cpus.lm2)$sigma)
```
