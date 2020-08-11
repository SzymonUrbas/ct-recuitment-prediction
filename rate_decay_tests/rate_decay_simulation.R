library(xtable)
SEED = 150
##### -------------------------------------------------------------
lambdas = c(5,10,20,50,100,200); R = c(1,0.9, 0.8, 0.7, 0.6, 0.5)
s = c(1,1); alpha = 0.05

## Likelihood-ratio test
set.seed(SEED)
power.tableLRT  = LRT.power.table(lambdas, s, R, alpha, iter = 5e6)

##### -------------------------------------------------------------
nums = 10 # number of intervals
lambdas = 1/nums*c(5,10,20,50,100,200); R = c(1,0.9, 0.8, 0.7, 0.6, 0.5)
s = c(nums, nums); alpha = 0.05

## Non-parametric bootstrapped test
set.seed(SEED)
power.tableNPBS  = NPBS.power.table(lambdas, s, R, alpha, iter = 5e4)



print(xtable(power.tableLRT));
print(xtable(power.tableNPBS))