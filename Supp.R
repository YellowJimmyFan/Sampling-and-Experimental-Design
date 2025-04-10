---
title: "Project Code"
output:
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r Defining preliminaries}
## Population data summaries
N1 <- 434; N2 <- 229; N3 <- 1439; N <- N1 + N2 + N3; Nhs <- c(N1, N2, N3)
W1 <- N1/N; W2 <- N2/N; W3 <- N3/N; Whs <- c(W1, W2, W3); n <- 325
## Function we use to find the total_estimate
Population_total_estimate <- function(N, Ws, ybars) {
  result <- N * sum(Ws*ybars)
  return(result)
}
## Function we use to find the proportion_estimate
Population_proportion_estimate <- function(Whs, ms, nhs) {
  Phs <- ms/nhs
  result <- sum(Whs*Phs)
  return(result)
}
## Function we use to find the variance of proportion under Proportional Alloc
Variance_proportional_prop <- function(N, n, Ws, sigmas) {
  result <- (1-n/N)*1/n*sum(Ws*sigmas)
  return(result)
}
## We can then find the variance of proportion under Neyman allocation by 
##   simply multiplying the Variance_proportional_prop result by N^2.
## Function we use to find the variance of proportion under Neyman allocation
Variance_neyman_prop <- function(n, Ws, sigmas, N) {
  result <- 1/n*(sum(Ws*sqrt(sigmas)))^2-1/N*sum(Ws*sigmas)
  return(result)
}
## We can then find the variance of total under Neyman allocation by simply
##   multiplying the Variance_neyman_prop result by N^2.
```

```{r Application under proportional allocation}
n1_prop <- 67; n2_prop <- 35; n3_prop <- 223
nhs_prop <- c(n1_prop, n2_prop, n3_prop)
proportional_samples <- read.csv("Proportional Allocation Data.csv")
Stats.Med_sample_prop <- proportional_samples[1:67, ]
Med.Res_sample_prop <- proportional_samples[68:102, ]
J.Clin_sample_prop <- proportional_samples[103:325, ]
## Calculating ybars
y1bar.prop <- mean(Stats.Med_sample_prop$Meta.Analysis)
y2bar.prop <- mean(Med.Res_sample_prop$Meta.Analysis)
y3bar.prop <- mean(J.Clin_sample_prop$Meta.Analysis)
ybars.prop <- c(y1bar.prop, y2bar.prop, y3bar.prop)
ybars_table_prop <- data.frame(
  y1bar.prop = y1bar.prop,
  y2bar.prop = y2bar.prop,
  y3bar.prop = y3bar.prop
); ybars_table_prop
## Sample errors
s1.prop <- var(Stats.Med_sample_prop$Meta.Analysis)
s2.prop <- var(Med.Res_sample_prop$Meta.Analysis)
s3.prop <- var(J.Clin_sample_prop$Meta.Analysis)
sample_errors.prop <- c(s1.prop, s2.prop, s3.prop)
samplerror_table_prop <- data.frame(
  s1.prop = s1.prop, s2.prop = s2.prop, s3.prop = s3.prop
); samplerror_table_prop
## Sample totals mi
m1.prop <- sum(Stats.Med_sample_prop$Meta.Analysis)
m2.prop <- sum(Med.Res_sample_prop$Meta.Analysis)
m3.prop <- sum(J.Clin_sample_prop$Meta.Analysis)
mhs.prop <- c(m1.prop, m2.prop, m3.prop)
sampltotal_table_prop <- data.frame(
  m1 = m1.prop, m2 = m2.prop, m3 = m3.prop
); sampltotal_table_prop
## Calculating the results
prop_total.est <- Population_total_estimate(N, Whs, ybars.prop)
prop_prop.est <- Population_proportion_estimate(Whs, mhs.prop, nhs_prop)
prop_prop.var <- Variance_proportional_prop(N, n, Whs, sample_errors.prop)
prop_total.var <- N^2*prop_prop.var
prop_total.se <- sqrt(prop_total.var); prop_prop.se <- sqrt(prop_prop.var)
## Results summary and confidence intervals
total_summary_prop <- data.frame(
  total_estimated_value = c(prop_total.est, NA),
  total_variance = c(prop_total.var, NA),
  total_standard_error = c(prop_total.se, NA),
  Confidence_Interval = prop_total.est + c(-1,1)*qnorm(0.975)*prop_total.se
); total_summary_prop
proportional_summary_prop <- data.frame(
  proportion_estimated_value = c(prop_prop.est, NA),
  proportion_variance = c(prop_prop.var, NA),
  proportion_standard_error = c(prop_prop.se, NA),
  Confidence_Interval = prop_prop.est + c(-1,1)*qnorm(0.975)*prop_prop.se
); proportional_summary_prop
```
```{r warning=FALSE}
## Strata summaries, Population variances
Stats.Med <- read.csv("Stats. Med. Data.csv")
Med.Res <- read.csv("Stat. Methods Med. Res Data.csv")
J.Clin <- read.csv("J. Clin. Epidemiol Data.csv")
## Data summaries of Stats.Med
Stats.Med_mean <- mean(Stats.Med$Meta.Analysis)
Stats.Med_total <- sum(Stats.Med$Meta.Analysis)
Stats.Med_variance <- var(Stats.Med$Meta.Analysis)
Stats.Med_sd <- sqrt(Stats.Med_variance)
## Data summaries of Med.Res
Med.Res_mean <- mean(Med.Res$Indicator); Med.Res_total <- sum(Med.Res$Indicator)
Med.Res_variance <- var(Med.Res$Indicator); Med.Res_sd <- sqrt(Med.Res_variance)
## Data summaries of J.Clin
J.Clin_mean <- mean(J.Clin$Indicator); J.Clin_total <- sum(J.Clin$Indicator)
J.Clin_variance <- var(J.Clin$Indicator); J.Clin_sd <- sqrt(J.Clin_variance)
## Overall summaries(order by Stats.Med, Med.Res, J.Clin)
strata_means <- c(Stats.Med_mean, Med.Res_mean, J.Clin_mean)
strata_totals <- c(Stats.Med_total, Med.Res_total, J.Clin_total)
strata_variances <- c(Stats.Med_variance, Med.Res_variance, J.Clin_variance)
strata_sds <- c(Stats.Med_sd, Med.Res_sd, J.Clin_sd)
## Show the results
stratameans_table <- data.frame(
  Stats.Med_mean = Stats.Med_mean, Med.Res_mean = Med.Res_mean,
  J.Clin_mean = J.Clin_mean
); stratameans_table
stratatotals_table <- data.frame(
  Stats.Med_total = Stats.Med_total, Med.Res_total = Med.Res_total,
  J.Clin_total = J.Clin_total
); stratatotals_table
stratavariances_table <- data.frame(
  Stats.Med_variance = Stats.Med_variance,
  Med.Res_variance = Med.Res_variance, J.Clin_variance = J.Clin_variance
); stratavariances_table
stratasd_table <- data.frame(
  Stats.Med_sd = Stats.Med_sd, Med.Res_sd = Med.Res_sd, J.Clin_sd = J.Clin_sd
); stratasd_table
```
```{r Neyman allocation sample size calculation}
## Calculate the size for neyman allocation
neyman_size_calc <- function(n, Ws, sigmas) {
  result <- n*Ws*sqrt(sigmas)/sum(Ws*sqrt(sigmas))
  return(result)
}
neyman_sizes <- neyman_size_calc(n, Whs, strata_variances); neyman_sizes
```
```{r Application of function under Neyman Allocation}
n1_neyman <- 97; n2_neyman <- 38; n3_neyman <- 190
nhs_neyman <- c(n1_neyman, n2_neyman, n3_neyman)
neyman_samples <- read.csv("Neyman Allocation Data.csv")
## Journals as Strata
Stats.Med_sample_neyman <- neyman_samples[1:97, ]
Med.Res_sample_neyman <- neyman_samples[98:135, ]
J.Clin_sample_neyman <- neyman_samples[136:325, ]
## Calculating ybars
y1bar.neyman <- mean(Stats.Med_sample_neyman$Meta.Analysis)
y2bar.neyman <- mean(Med.Res_sample_neyman$Meta.Analysis)
y3bar.neyman <- mean(J.Clin_sample_neyman$Meta.Analysis)
ybars.neyman <- c(y1bar.neyman, y2bar.neyman, y3bar.neyman)
ybars_table_neyman <- data.frame(
  y1bar.neyman = y1bar.neyman, y2bar.neyman = y2bar.neyman, 
  y3bar.neyman = y3bar.neyman
); ybars_table_neyman
## Sample errors
s1.neyman <- var(Stats.Med_sample_neyman$Meta.Analysis)
s2.neyman <- var(Med.Res_sample_neyman$Meta.Analysis)
s3.neyman <- var(J.Clin_sample_neyman$Meta.Analysis)
sample_errors.neyman <- c(s1.neyman, s2.neyman, s3.neyman)
samplerror_table_neyman <- data.frame(
  s1.neyman = s1.neyman, s2.neyman = s2.neyman, s3.neyman = s3.neyman
); samplerror_table_neyman
## Sample totals mi
m1.neyman <- sum(Stats.Med_sample_neyman$Meta.Analysis)
m2.neyman <- sum(Med.Res_sample_neyman$Meta.Analysis)
m3.neyman <- sum(J.Clin_sample_neyman$Meta.Analysis)
mhs.neyman <- c(m1.neyman, m2.neyman, m3.neyman)
sampltotal_table_neyman <- data.frame(
  m1 = m1.neyman, m2 = m2.neyman, m3 = m3.neyman
); sampltotal_table_neyman
## Calculating the results
neyman_total.est <- Population_total_estimate(N, Whs, ybars.neyman)
neyman_prop.est <- Population_proportion_estimate(Whs, mhs.neyman, nhs_neyman)
neyman_prop.var <- Variance_neyman_prop(n, Whs, strata_variances, N)
neyman_total.var <- N^2*neyman_prop.var; neyman_total.se <- sqrt(neyman_total.var)
neyman_prop.se <- sqrt(neyman_prop.var)
## Results summary and confidence intervals
total_summary_neyman <- data.frame(
  total_estimated_value = c(neyman_total.est, NA),
  total_variance = c(neyman_total.var, NA),
  total_standard_error = c(neyman_total.se, NA),
  Confidence_Interval = neyman_total.est + c(-1,1)*qnorm(0.975)*neyman_total.se
); total_summary_neyman
neyman_summary_prop <- data.frame(
  proportion_estimated_value = c(neyman_prop.est, NA),
  proportion_variance = c(neyman_prop.var, NA),
  proportion_standard_error = c(neyman_prop.se, NA),
  Confidence_Interval = neyman_prop.est + c(-1,1)*qnorm(0.975)*neyman_prop.se
); neyman_summary_prop
```

