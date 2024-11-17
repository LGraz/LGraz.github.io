---
title: |
  Linear Mixed Models
  Simplifications Tips and Tricks
author: Lukas Graz
note-type: quarto
date: 2024-07-26T00:00:00.000Z
uuid: de176b8c-de72-4a9b-8705-3ccb144be6cf
permalink: /posts/simplify-lmm-ticks-tricks
execute:
  cache: false
  message: false
  warning: false
  toc: true
tags:
  - simulation
  - power
  - lmm
---


1.  ~~Technical Details~~
2.  Sleepstudy
3.  Finger Tapping
4.  Power Simulations

# Sleepstudy

**Reaction ~ Days**

. . .

``` r
library(lmerTest)
library(simr)

data(sleepstudy)
S <- sleepstudy
```

![](../quarto/lmer-slopes-new/lmer_vs_ttests-new_files/figure-commonmark/plot-sleepstudy-1.png)

## Model Illustration

![source: bookdown.org/steve_midway/DAR](../assets/images/07_models.png)

## Model Equation

``` r
fit <- lmer(Reaction ~ Days + (Days|Subject), S)
```

. . .

``` r
equatiomatic::extract_eq(fit)
```

$$
\begin{aligned}
  \operatorname{Reaction}_{i}  &\sim N \left(\alpha_{j[i]} + \beta_{1j[i]}(\operatorname{Days}), \sigma^2 \right) \\    
\left(
  \begin{array}{c} 
    \begin{aligned}
      &\alpha_{j} \\
      &\beta_{1j}
    \end{aligned}
  \end{array}
\right)
  &\sim N \left(
\left(
  \begin{array}{c} 
    \begin{aligned}
      &\mu_{\alpha_{j}} \\
      &\mu_{\beta_{1j}}
    \end{aligned}
  \end{array}
\right)
, 
\left(
  \begin{array}{cc}
     \sigma^2_{\alpha_{j}} & \rho_{\alpha_{j}\beta_{1j}} \\ 
     \rho_{\beta_{1j}\alpha_{j}} & \sigma^2_{\beta_{1j}}
  \end{array}
\right)
 \right)
    \text{, for Subject j = 1,} \dots \text{,J}
\end{aligned}
$$

## Test

**LMM**

``` r
coef(summary(fit))["Days", , drop=F]
```

    #>      Estimate Std. Error df t value  Pr(>|t|)
    #> Days    10.47      1.546 17   6.771 3.264e-06

. . .

**Aggregate + t-Test**

``` r
COEFS <- coef(lmList(Reaction ~ Days | Subject, S))
head(COEFS)
```

    #>     (Intercept)   Days
    #> 308       244.2 21.765
    #> 309       205.1  2.262
    #> 310       203.5  6.115
    #> 330       289.7  3.008
    #> 331       285.7  5.266
    #> 332       264.3  9.567

. . .

``` r
t.test(COEFS$Days)$p.value
```

    #> [1] 3.264e-06

# Finger Tapping

**Tap_Speed ~ Fatigued + Caffein + Physical_Activity + Time**

. . .

![Consultancy Case: Sarah Stadelmann
(Illustration)](../assets/images/example-illustration.png)

## LMM

``` r
D <- readRDS("./FingerTapping.rds")
anova(RSLP_fit <- 
  lmer(Tap_Speed_t0 ~ 
     Fatigued + Caffein + Physical_Activity + Time +
    (Fatigued + Caffein + Physical_Activity + Time || PersonID), D
))
```

    #> Type III Analysis of Variance Table with Satterthwaite's method
    #>                   Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)    
    #> Fatigued           1.118   1.118     1  20.3   16.52 0.00059 ***
    #> Caffein            0.710   0.710     1  16.7   10.48 0.00492 ** 
    #> Physical_Activity  0.134   0.134     1  14.3    1.97 0.18153    
    #> Time               0.843   0.843     1  27.5   12.45 0.00149 ** 
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Aggregate + t-Test

1.  `lm` per person
2.  Extract coefficients
3.  t-test

. . .

``` r
COEFS <- coef(lmList(Tap_Speed_t0 ~ 
  Fatigued + Caffein + Physical_Activity + Time| PersonID, D))
head(COEFS)
```

    #>   (Intercept)  Fatigued Caffein Physical_Activity      Time
    #> 1       5.537 -0.005006 0.11454           0.01294  0.010271
    #> 2       5.954  0.008884 0.06979           0.24502 -0.009993
    #> 4       6.089 -0.091132 0.24615                NA -0.006202
    #> 6       5.690 -0.002435 0.09211                NA -0.025624
    #> 7       6.197 -0.139987      NA                NA  0.034949
    #> 8       6.288 -0.042551      NA           0.45775  0.010620

. . .

``` r
t.test(COEFS$Fatigued)$p.value
```

    #> [1] 0.006832

## Comparison

<div class="columns">

<div class="column" width="50%">

**LMM**

- Challenging numerics  
- Crossed random effects  

</div>

<div class="column" width="50%">

**aggregate + t-Test**

- KISS  
- Fewer assumptions  

</div>

</div>

## Type I Error & Power

``` r
##| code-fold: true
##| code-summary: "Simulation Code"

# Setup auxiliary functions
get_type1_CI <- function(P.vals){
  C.int <- apply(P.vals < 0.05, 1, function(x)
    prop.test(sum(x), length(x), p = 0.05)$conf.int |> round(4))

  cbind(
    `Power` = rowMeans(P.vals < 0.05),
    lwr = C.int[1, ],
    upr = C.int[2, ]
  ) |> as.data.frame()
}

get_p_vals <- function(Y){
  D$Y <- Y
  c(
    LMM = lmer(Y ~ 
       Fatigued + Caffein + Physical_Activity + Time +
      (Fatigued + Caffein + Physical_Activity + Time || PersonID), D) |> 
      summary() |> coef() |> _["Fatigued", "Pr(>|t|)"],
    `t-Test` = lmList(Y ~ 
      Fatigued + Caffein + Physical_Activity + Time | PersonID, D
      ) |> coef() |> _$Fatigued |> t.test() |> _$p.value
  )
}

# Simulate H0 ===========
xfun::cache_rds({
    RSLP_fit <- lmer(Tap_Speed_t0 ~ 
        Fatigued + Caffein + Physical_Activity + Time +
       (Fatigued + Caffein + Physical_Activity + Time | PersonID), D
    )
    fixef(RSLP_fit) <- c(0,0,0,0,0)
    simulate(
      RSLP_fit, 
      nsim = 10000, 
      seed=123) |> 
      parallel::mclapply(get_p_vals) 
  },
  file = "cache_sim3.1.rds",
  dir = "cache/",
  hash = list(D, as.character(body(get_p_vals)))
) |> simplify2array() |> get_type1_CI() -> H0_CI

# Simulate H1 ===========
xfun::cache_rds({
    RSLP_fit <- lmer(Tap_Speed_t0 ~ 
        Fatigued + Caffein + Physical_Activity + Time +
       (Fatigued + Caffein + Physical_Activity + Time | PersonID), D
    )
    fixef(RSLP_fit) <- c(0,0.02,0,0,0)
    simulate(
      RSLP_fit, 
      nsim = 10000, 
      seed=123) |> 
      parallel::mclapply(get_p_vals) 
  },
  file = "cache_sim3.2.rds",
  dir = "cache/",
  hash = list(D, as.character(body(get_p_vals)))
) |> simplify2array() |> get_type1_CI() -> H1_CI

# Plot results ============
H0_CI$Hypothesis <- "H0: Fatigue = 0"
H1_CI$Hypothesis <- "H1: Fatigue = 0.02"

H_all <- rbind(H0_CI, H1_CI) 
H_all$test <- rownames(H_all)

ggplot(H_all, aes(x = test, y = `Power`, ymin = lwr, ymax = upr)) +
  geom_pointrange() +
  geom_hline(yintercept = 0.05, linetype = 2) +
  facet_wrap(~Hypothesis, scales = "free") +
  theme_minimal() +
  xlab("")
```

![](../quarto/lmer-slopes-new/lmer_vs_ttests-new_files/figure-commonmark/Ex-sim-1.png)

# Power-Simulations

## Using `simr`

``` r
powerSim(RSLP_fit)
```

. . .

``` r
powerCurve(RSLP_fit, along = "PersonID") |> plot()
```

![](../quarto/lmer-slopes-new/lmer_vs_ttests-new_files/figure-commonmark/powercurve_sim-1.png)

## Simulate with `rnorm()` & co.

![](../assets/images/boring.png)

## `lme4::simulate.merMod`

1.  Existing model
2.  Modify existing model (`simr`)

. . .

- `fixef<-`
- `coef<-`
- `VarCorr<-`
- `sigma<-`
- `scale<-`

. . .

``` r
fixef(RSLP_fit) <- c(0,0,0,0,0)
simulate(RSLP_fit, nsim = 100)
```

## No existing model?

`simulate(Y ~ ..., newparams = list(theta =  ???, ...))`[^1]

. . .

``` r
simr::makeLmer(
  Y ~ Fatigued + Caffein + Physical_Activity + Time +
     (Fatigued + Caffein + Physical_Activity + Time | PersonID), 
  fixef = c(0,0,0,0,0),
  VarCorr = diag(5),
  data = D,
  sigma = 1
)
```

    #> Linear mixed model fit by REML ['lmerMod']
    #> Formula: Y ~ Fatigued + Caffein + Physical_Activity + Time + (Fatigued +  
    #>     Caffein + Physical_Activity + Time | PersonID)
    #>    Data: D
    #> REML criterion at convergence: 2018
    #> Random effects:
    #>  Groups   Name              Std.Dev. Corr               
    #>  PersonID (Intercept)       1                           
    #>           Fatigued          1        0.00               
    #>           Caffein           1        0.00 0.00          
    #>           Physical_Activity 1        0.00 0.00 0.00     
    #>           Time              1        0.00 0.00 0.00 0.00
    #>  Residual                   1                           
    #> Number of obs: 567, groups:  PersonID, 27
    #> Fixed Effects:
    #>       (Intercept)           Fatigued            Caffein  Physical_Activity  
    #>                 0                  0                  0                  0  
    #>              Time  
    #>                 0

## Caching: `xfun::cache_rds`

``` r
xfun::cache_rds({
    expression
  },
  hash = list_with_relevant_objects
)
```

## Have Fun Simulating!

1.  Aggregate per Subject $\rightarrow$ t-Test
2.  `equatiomatic::extract_eq()`
3.  `lme4::simulate.merMod()`
    1.  `fixef<-` & co.
    2.  `simr::makeLmer()`
4.  `xfun::cache_rds()`

[^1]: theta: cholesky factor of the normalized random effects covariance
    matrix.
