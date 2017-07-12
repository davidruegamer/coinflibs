# coinflibs
**Co**nditional **Inf**erence after **Li**kelihood-**b**ased **S**election

The package contains functions to calculate limits and conduct inference in a selective manner for linear models after likelihood- or test-based model selection. 

Please see https://arxiv.org/abs/1706.09796 for more details.

Future improvements:
  - Testing grouped baselearners
  - Correct for plugin error variance approach
 
 
## Example
 
```R
# install and load package
install_github("davidruegamer/coinflibs")
library("coinflibs")

library(MASS)
# Fit initial model
cpus$perf <- log10(cpus$perf)
mod <- lm(perf ~ .-name, data = cpus)

# use the stepAIC function to find the best model in a backward 
# stepwise search
cpus.lm <- stepAIC(mod, trace = FALSE, direction = "backward")

# recalculate all visited models
lom <- c(lapply(colnames(model.matrix(mod)), function(x) 
  update(mod, as.formula(paste0("perf ~ .-", x)))), list(mod))
  
# extract the components of all visited models
compsAIC <- extract_components(listOfModels = lom,
                               response = cpus$perf,
                               what = c("aic"))
                                 
# perform likelihood ratio test at level 
alpha = 0.001

# check for non-significant variables
coefTable <- summary(cpus.lm)$coefficients
drop <- rownames(coefTable)[alpha < coefTable[,4]]

# drop non-significant variable
cpus.lm2 <- update(cpus.lm, as.formula(paste0(".~.-",drop)))

# extract components associated with the LRT comparison
compsLRT <- extract_components(list(cpus.lm, cpus.lm2),
                               response = cpus$perf, 
                               what = "lrt",
                               alpha = alpha)
                               
# naive inference
unadj_pvs <- summary(cpus.lm2)$coefficients[,4]

# now extract testvector, calculate limits and perform selective tests
# test vectors

vTs <- extract_testvec(cpus.lm2)

# calculate limits
limitsAIC <- calculate_limits(compsAIC, vTs)
limitsLRT <- calculate_limits(compsLRT, vTs)

# check restriction on p-values separately
cbind(
selinf(limitObject = limitsAIC, y = cpus$perf, sd = sigma(cpus.lm2)),
unadjusted_pval = unadj_pvs
)

cbind(
selinf(limitObject = limitsLRT, y = cpus$perf, sd = sigma(cpus.lm2)),
unadjusted_pval = unadj_pvs
)

# calculate p-values (does automatically combine limitObjects)
res <- selinf(limitObject = list(limitsAIC, limitsLRT), 
                       y = cpus$perf, 
                       # plugin estimate for true value
                       sd = sigma(cpus.lm2))
                       
cbind(res, unadjusted_pval = unadj_pvs)
```
