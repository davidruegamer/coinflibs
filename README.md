# coinflibs
**Co**nditional **Inf**erence after **Li**kelihood-**b**ased **S**election

The package contains functions to calculate limits and conduct inference in a selective manner for linear models after likelihood- or test-based model selection. 

Please see https://arxiv.org/abs/1706.09796 for more details.

## Example: Combining stepwise AIC search and significance hunting
 
```R
# install and load package
install_github("davidruegamer/coinflibs")
library("coinflibs")

library(MASS)
# use the cpus data
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
allvisited <- lapply(attr(mod$terms, "term.labels"), 
                     function(x) update(mod, as.formula(paste0("perf ~ .-", x))))
# combine the visited models and the final model                     
lom1 <- c(allvisited, list(mod))

# perform F-test at level 
alpha = 0.001
# and check for non-significant variables
coefTable <- anova(cpus.lm)
drop <- rownames(coefTable)[alpha < coefTable[-nrow(coefTable),5]]

# drop non-significant variable
cpus.lm2 <- update(cpus.lm, as.formula(paste0(".~.-",drop)))

# create list of models, which are examined during the significance search
lom2 <- list(cpus.lm, cpus.lm2)

# now compute selective inference, which adjust for the AIC-based search as
# well as significance hunting
selinf( # supply all lists of visited models, where the best model in the 
        # first list is interpreted as the final model
       lom2, # list given by signficance hunting
       lom1, # list given by AIC-based search
       response = cpus$perf,
       what = c("Ftest", "aic"), # specify what type of selection was done 
                                 # for each supplied list
       sd = summary(cpus.lm2)$sigma
      )
```
