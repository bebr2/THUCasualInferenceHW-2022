

```{r}
getwd()
setwd("D:/个人/学习/因果推断导论/大作业")
```
```{r}
data <- read.csv("./data/mydata.csv")
head(data)
attach(data)
barplot(data[1:51,]$CIPERTHOUSAND)
data$SN[which.max(data$CIPERTHOUSAND)]
lbs1 = c("已经废死", "未废死")
co = c(sum(data$ABOLISHDP[data$ABOLISHDP==1]),102-sum(data$ABOLISHDP[data$ABOLISHDP==1]))
lbs2 = paste(lbs1, " ", co, "个州")
pie(c(sum(data$ABOLISHDP[data$ABOLISHDP==1]), 102-sum(data$ABOLISHDP[data$ABOLISHDP==1])), labels=lbs2)
```
```{r}
fit1 <- lm(CIPERTHOUSAND~ABOLISHDP, data=data)
summary(fit1)
plot(ABOLISHDP, CIPERTHOUSAND)
abline(fit1)
```
```{r}
x0 <- as.matrix(data[data$ABOLISHDP==0,3:13])
x1 <- as.matrix(data[data$ABOLISHDP==1,3:13])
mahalanobis(colMeans(x1), colMeans(x0), (var(x0)+var(x1))/2)

```
```{r}
plot_difference = function(i, name)  {
  ymax <- max(c(density(x0[,i])$y,density(x1[,i])$y))
  plot(density(x0[,i]),ylim = c(0,1.1*ymax), main=name)
  lines(density(x1[,i]),col="red")
}
```


```{r}
plot_difference(2, "high_school_or_higher_of_over25")
plot_difference(4, "white")
plot_difference(9, "no_health_insurance")
```


```{r}
fit2 <- lm(CIPERTHOUSAND~.-X-SN-CIPERTHOUSAND,data=data)
summary(fit2)
```


```{r}
md1 <- glm(ABOLISHDP~NHIP+HSOHO25+UEO16+PP+BLACK, data = data, family = binomial)
summary(md1)
```


```{r}

fw_LRT <-function(object, scope, C) {
  step <- 0
  cov_name <- names(coef(scope))[-1]
  tab <- matrix(NA, nrow=length(cov_name),ncol=length(cov_name)- length(coef(object))+1)
  rownames(tab) <- cov_name
  repeat {
    if (setequal(names(coef(object)),names(coef(scope)))) {
      cat("\nWARNING: Maximun Model Reached.\n")
      break
    }
    step <- step+1
    cat("\n>>> Step =", step, "\n")
    res <-add1(object, scope = scope, test="LRT")
    print(res)
    tab[rownames(res)[-1], step] <- res$LRT[-1]
    add_id <-which.max(res$LRT)
    if(res$LRT[add_id]<C)
      break
    object <- update(object,as.formula(paste("~ . +",rownames(res)[add_id])))
    }
  cat("\n>>> Final Model\n")
  print(object)
  return(list(md = object, tab = tab[, 1:step]))
}

```


```{r}
mdL <- glm(ABOLISHDP~.-X-SN-CIPERTHOUSAND, data = data, family = binomial)
res <-fw_LRT(md1, scope = mdL, C = 1)
```


```{r}
summary(res$md)
```


```{r}
print(round(res$tab, 1), na.print = "-")
```


```{r}
md2 <- res$md
summary(md2)
```


```{r}
term <-names(coef(md2))[-1]
mdQ <-glm(as.formula(paste("ABOLISHDP ~", "(",paste(term, collapse = " + "),")^2 +",paste(paste0("I(", term, "^2)"), collapse = " + "))),data, family = binomial)
```


```{r}
res <-fw_LRT(md2, scope = mdQ, C=qnorm(0.975)^2)
```


```{r}
print(round(res$tab, 1), na.print = "-")
md3 <- res$md
summary(md3)
```


```{r}
plot(data$ABOLISHDP, predict(md2, type = "link"), xlab = "W", ylab = "lps",main = "Linearized Propensity Score")
```


```{r}
data$le = md3$linear.predictor
number = nrow(data)
number.t = sum(data$ABOLISHDP)
number.c = number - number.t
print(number.t)
print(number.c)
```


```{r}
lbound = min(data$le[data$ABOLISHDP==1])
ubound = max(data$le[data$ABOLISHDP==0])
print(lbound)
print(ubound)
data.trim = data[data$le >= lbound & data$le <= ubound, ]

number = nrow(data.trim)
number.t = sum(data.trim$ABOLISHDP)
number.c = number - number.t
tstat = (mean(data.trim$le[data.trim$ABOLISHDP==1]) - mean(data.trim$le[data.trim$ABOLISHDP==0]))/sqrt((var(data.trim$le[data.trim$ABOLISHDP==1])*(number.t - 1)
+ var(data.trim$le[data.trim$ABOLISHDP==0])*(number.c - 1))/(number - 2)*(1/number.c + 1/number.t))
cbind(lbound, ubound, number.c, number.t, tstat)
```


```{r}
data$le = md2$linear.predictor
lbound = min(data$le[data$ABOLISHDP==1])
ubound = max(data$le[data$ABOLISHDP==0])
print(lbound)
print(ubound)
data.trim = data[data$le >= lbound & data$le <= ubound, ]

number = nrow(data.trim)
number.t = sum(data.trim$ABOLISHDP)
number.c = number - number.t
tstat = (mean(data.trim$le[data.trim$ABOLISHDP==1]) - mean(data.trim$le[data.trim$ABOLISHDP==0]))/sqrt((var(data.trim$le[data.trim$ABOLISHDP==1])*(number.t - 1)
+ var(data.trim$le[data.trim$ABOLISHDP==0])*(number.c - 1))/(number - 2)*(1/number.c + 1/number.t))
cbind(lbound, ubound, number.c, number.t, tstat)
```


```{r}
data.trim$b4 = 1
library(Hmisc)
data.trim$b4[data.trim$ABOLISHDP==1] = as.numeric(cut2(data.trim$le[data.trim$ABOLISHDP==1], g = 4))
le.t = sort(data.trim$le[data.trim$ABOLISHDP==1])
data.trim$b4[data.trim$ABOLISHDP==0] = cut(data.trim$le[data.trim$ABOLISHDP==0], breaks = c(-100, le.t[8], le.t[16], le.t[24], 100))
sum(data.trim$b4==1)-8
sum(data.trim$b4==2)-8
sum(data.trim$b4==3)-8
sum(data.trim$b4==4)-8
```


```{r}
balan = NULL
for(j in 1:4){
  for(i in 0:1) balan = cbind(balan, apply(data.trim[data.trim$ABOLISHDP==i & data.trim$b4==j, c(3:13, 16)], 2, mean))
}
colnames(balan) = c("1c", "1t", "2c", "2t", "3c", "3t", "4c", "4t")
round(balan, 2)
```



```{r}
lbound = min(data$le[data$ABOLISHDP==1])
ubound = max(data$le[data$ABOLISHDP==0])
print(lbound)
print(ubound)
data2 = data[data$le >= lbound & data$le <= ubound, ]
library(Matching)

match_ps <-Match(Y = data2$CIPERTHOUSAND, Tr = data2$ABOLISHDP, X = data2$le, M = 1,replace = FALSE, distance.tolerance = 0)
summary(match_ps)
```


```{r}
MatchBalance(data2$ABOLISHDP ~ ., data = data2[,3:13], match.out = match_ps)
```


```{r}
get_table <- function(X, W, lps, basic = TRUE, alpha = 0.05,
left.include = FALSE) {
  K <- ncol(X)
  ps <- exp(lps)/(1 + exp(lps)) # propensity score in [0, 1]
  ind <- list("0" = which(W == 0), "1" = which(W == 1))
  X_aug <- cbind(X, NA, ps, lps)
  colnames(X_aug)[K + (1:3)] <- c("Multivar_measure", "pscore", "linear_pscore")
  tab <- t(apply(X_aug, 2, function(x) # basic summary statistics
  unlist(tapply(x, W, function(z) c(mean(z), sd(z))))))
  mu_dif <- tab[, 3] - tab[, 1] # mean(Xt) - mean(Xc)
  NorDif <- mu_dif/sqrt((tab[, 2]^2 + tab[, 4]^2)/2) # normalized difference
  LR_STD <- log(tab[, 4]) - log(tab[, 2]) # Log ratio of standard deviation
  get_pi <- function(x, grp) { # measure M3 & M4
  if (anyNA(x))
  return(NA)
  else {
  Fun <- function(xx, boundary = TRUE) {
  if (boundary)
  return(mean(x[ind[[grp + 1]]] <= xx))
  else
  return(mean(x[ind[[grp + 1]]] < xx))
  }
  qtl <- quantile(x[ind[[2 - grp]]], c(1 - alpha/2, alpha/2), type = 1)
  return(1 - Fun(qtl[1]) + Fun(qtl[2], boundary = left.include))
  }
  }
  pi_c <- apply(X_aug, 2, get_pi, grp = 0)
  pi_t <- apply(X_aug, 2, get_pi, grp = 1)
  # Mahalanobis distance
  Sigma_ct <- (cov(X[ind$"0", ]) + cov(X[ind$"1", ]))/2
  mv_ct <- sqrt(t(mu_dif[1:K]) %*% (Sigma_ct) %*% mu_dif[1:K])
  res <- cbind(NorDif, LR_STD, pi_c, pi_t) # merge tables
  if (basic) {
  colnames(tab) <- c("Mean_c", "(S.D.)_c", "Mean_t", "(S.D.)_t")
  res <- cbind(tab, res)
  }
  res[K + 1, ncol(res) - 3] <- mv_ct
  return(res)
}
```


```{r}
tab <- get_table(data2[, 3:13], data2$ABOLISHDP, data2$le)
```


```{r}
dat_ps <- data2[c(match_ps$index.treated, match_ps$index.control), ]
tab_ps <-get_table(dat_ps[,3:13], dat_ps$ABOLISHDP, data2$le[c(match_ps$index.treated,match_ps$index.control)], basic = FALSE)
tab_ps
```


```{r}
K <- ncol(data2) - 2
par(mfrow = c(1, 1))
plot(NULL, xlim = c(-2, 2), ylim = c(1, K), xlab = "Normalized Difference",
ylab = "Covariate Label")
abline(h = 1:K, lty = 2, col = "gray65")
abline(v = 0, col = "gray65")
points(tab[1:K, "NorDif"], 1:K, pch = 8)
points(tab_ps[1:K, "NorDif"], 1:K, pch = 1, col = "red")
```



```{r}
no_Covariate_tau = mean(data2[match_ps$index.treated,]$CIPERTHOUSAND) - mean(data2[match_ps$index.control,]$CIPERTHOUSAND)
no_Covariate_tau
no_Covariate_var = 2 * (var(data2[cbind(match_ps$index.treated, match_ps$index.control),]$CIPERTHOUSAND)) / (nrow(data2[match_ps$index.treated,]))
no_Covariate_var
```


```{r}
names(md2$coefficients[-1])
```




```{r}
library(MASS)
match_data = data2[cbind(match_ps$index.treated, match_ps$index.control),3:15]
lmo<-lm(CIPERTHOUSAND~ABOLISHDP+UEO16+PP,data=match_data)
lmfor<-step(lmo,scope=list(upper=~ABOLISHDP+WHITE+BLACK+PP+UEO16+HSOHO25+NHIP+LT1WH+HSOH1824+BOHO25+PF+MIH,lower=~ABOLISHDP+UEO16+PP),direction="both")
summary(lmfor)
```


```{r}
stj1 = c(0.147101, 0.574112, 0.049618, 0.039387, 0.172255, 0.731092, 0.099568)
```




```{r}
Covariate_tau = 0.004128
Covariate_tau
Covariate_var = sqrt(sum(stj1^2)/7)
Covariate_var
```


```{r}
match_ps2 <-Match(Y = data2$CIPERTHOUSAND, Tr = data2$ABOLISHDP, X = data2$le, M = 1,replace = TRUE, distance.tolerance = 0)
tab2 <- get_table(data2[, 3:13], data2$ABOLISHDP, data2$le)
dat_ps2 <- data2[c(match_ps2$index.treated, match_ps2$index.control), ]
tab_ps2 <-get_table(dat_ps2[,3:13], dat_ps2$ABOLISHDP, data2$le[c(match_ps2$index.treated,match_ps2$index.control)], basic = FALSE)
tab_ps2
match_ps2$index.treated
match_ps2$index.control
```


```{r}
K <- ncol(data2) - 2
par(mfrow = c(1, 1))
plot(NULL, xlim = c(-2, 2), ylim = c(1, K), xlab = "Normalized Difference",
ylab = "Covariate Label")
abline(h = 1:K, lty = 2, col = "gray65")
abline(v = 0, col = "gray65")
points(tab2[1:K, "NorDif"], 1:K, pch = 8)
points(tab_ps[1:K, "NorDif"], 1:K, pch = 1, col = "red")
points(tab_ps2[1:K, "NorDif"], 1:K, pch = 3, col = "blue")
```


```{r}
data2[match_ps2$index.treated,]$CIPERTHOUSAND
data2[match_ps2$index.control,]$CIPERTHOUSAND
```


```{r}
no_Covariate_tau2 = mean(data2[match_ps2$index.treated,]$CIPERTHOUSAND) - mean(data2[match_ps2$index.control,]$CIPERTHOUSAND)
no_Covariate_tau2
all_var = var(data2[cbind(match_ps2$index.treated, match_ps2$index.control),]$CIPERTHOUSAND)
indicator<-duplicated(match_ps2$index.control)
data2[match_ps2$index.control,]$SN
table(match_ps2$index.control[indicator])+1
```


```{r}
no_duplicated = 32 + 32 - 23
no_Covariate_var2 = all_var * (no_duplicated +2^2+2^2+4^2+2^2+2^2+11^2)/ (32^2)
no_Covariate_var2


match_data2 = data2[cbind(match_ps2$index.treated, match_ps2$index.control),3:15]
lmo2<-lm(CIPERTHOUSAND~ABOLISHDP+UEO16+PP,data=match_data2)
lmfor2<-step(lmo2,scope=list(upper=~ABOLISHDP+WHITE+BLACK+PP+UEO16+HSOHO25+NHIP+LT1WH+HSOH1824+BOHO25+PF+MIH,lower=~ABOLISHDP+UEO16+PP),direction="both")
summary(lmfor2)
```


```{r}
lmfor2$coefficients[3:9]
stj2 = c(0.108004, 0.106363, 0.026918,0.097420,0.109739,0.049959,0.134507)
D_bar2 = colMeans(data2[match_ps2$index.treated,c(8,13,7,11,4,5,3)]) - colMeans(data2[match_ps2$index.control,c(8,13,7,11,4,5,3)])
D_bar2
stj2
```


```{r}
Covariate_tau2 = no_Covariate_tau2 - D_bar2 %*% c(0.2614760, 0.5910031, 0.2978525, 0.4287860, 0.3776103, 0.1150685, 0.1995873 )
Covariate_tau2
Covariate_var2 = sqrt(sum(stj2^2))/7
Covariate_var2
```





