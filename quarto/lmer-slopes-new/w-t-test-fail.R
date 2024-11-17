# Functions ======================================
w.t.test <- function(x, weights, data = NULL){
  lm(x ~ 1, data = data, weights = weights) |> 
    summary() |> coef() |> _["(Intercept)", "Pr(>|t|)"]
}

get_type1_CI <- function(P.vals){
  C.int <- apply(P.vals < 0.05, 1, function(x)
    prop.test(sum(x), length(x), p = 0.05)$conf.int |> round(3))

  rbind(
    `type I error` = rowMeans(P.vals < 0.05),
    lwr = C.int[1, ],
    upr = C.int[2, ]
  ) |> knitr::kable()
}

# Analysis ======================================
library(lmerTest)
library(simr)

# Subset
# only keep days 2:4 for the first 10 subjects
S.sub <- subset(sleepstudy, as.numeric(Subject) > 10 | Days %in% 2:4)
subj_n <- table(S.sub$Subject); subj_n

get.pvals <- function(y){
  # calculate the Subject intercepts
  COEFS <- lm(y ~ 0 + Subject, S.sub) |> 
    summary() |> coef() |> as.data.frame()
  itcpts    <- COEFS$Estimate
  itcpts.se <- COEFS$`Std. Error`

  # return p-values via different methods
  return(c(
    lmer          = lmer(y ~ 1 + (1|Subject), S.sub) |> 
                     summary() |> coef() |> _["(Intercept)", "Pr(>|t|)"],
    t.test         = w.t.test(itcpts, NULL),
    w.t.test.var   = w.t.test(itcpts, 1 / itcpts.se^2),
    w.t.test.var2  = w.t.test(itcpts, 1 /(itcpts.se^2 + var(itcpts))),
    w.t.test.subj_n= w.t.test(itcpts, subj_n)
  ))
}

get.pvals(S.sub$Reaction)

# Simulate H0 =================================
fit <- lmer(Reaction ~ 1 + (1|Subject), S.sub)
fixef(fit) <- 0  
sigma(fit)
sigma(fit) <- 100
fit@theta
simulate(fit, nsim = 10000, seed=123) |> 
  sapply(get.pvals) |> get_type1_CI()
#> |             |   lmer| t.test| w.t.test.var| w.t.test.var2| w.t.test.subj_n|
#> |:------------|------:|------:|------------:|-------------:|---------------:|
#> |type I error | 0.0438| 0.0478|        0.055|        0.0437|           0.055|
#> |lwr          | 0.0400| 0.0440|        0.051|        0.0400|           0.051|
#> |upr          | 0.0480| 0.0520|        0.060|        0.0480|           0.060|
#=>> weighted t.test not working as expected
