## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(hdcuremodels)
library(survival)

## -----------------------------------------------------------------------------
data <- generate_cure_data(N = 200, J = 50, nTrue = 10, A = 1.8, rho = 0.2)
training <- data$Training
testing <- data$Testing

## -----------------------------------------------------------------------------
km.train <- survfit(Surv(cryr, relapse.death) ~ 1, data = amltrain)

## ----echo=FALSE, fig=TRUE-----------------------------------------------------
plot(km.train, mark.time = TRUE, xlab = "Time (years)", ylab = "Relapse-free survival")

## ----echo=TRUE----------------------------------------------------------------
nonzerocure_test(km.train)

## -----------------------------------------------------------------------------
cure_estimate(km.train)

## ----echo=TRUE----------------------------------------------------------------
sufficient_fu_test(km.train)

## ----args---------------------------------------------------------------------
args(curegmifs)

## ----eval=FALSE---------------------------------------------------------------
#  fitgmifs <- curegmifs(Surv(cryr, relapse.death) ~ ., data = amltrain,
#                        x.latency = amltrain, model = "weibull")

## ----args2--------------------------------------------------------------------
args(cureem)

## -----------------------------------------------------------------------------
fitem <- cureem(Surv(cryr, relapse.death) ~ ., data = amltrain, 
                      x.latency = amltrain, model = "cox", 
                      lambda.inc=0.009993, lambda.lat=0.02655)

## ----args3--------------------------------------------------------------------
args(cv_cureem)

## -----------------------------------------------------------------------------
fit.cv <- cv_cureem(Surv(Time, Censor) ~ ., data = training,
                     x.latency = training, fdr.control = FALSE,
                     grid.tuning = FALSE, nlambda.inc = 10, nlambda.lat = 10,
                     n_folds = 2, seed = 23, verbose = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  lambda.inc <- lambda.lat <- rep(0, 100)
#  for (k in 1:100) {
#    print(k)
#    coxem.auc.k<-cv_cureem(Surv(cryr, relapse.death) ~ ., data = amltrain,
#                  x.latency = amltrain, model = "cox", penalty = "lasso",
#                  scale = TRUE, grid.tuning = TRUE, nfolds = 10, nlambda.inc = 20,
#                  nlambda.lat = 20, verbose = FALSE, parallel = TRUE, measure.inc = "auc")
#    lambda.inc[k]<-coxem.auc.k$selected.lambda.inc
#    coxem.c.k<-cv_cureem(Surv(cryr, relapse.death) ~ ., data = amltrain,
#                  x.latency = amltrain, model = "cox", penalty = "lasso",
#                  scale = TRUE, grid.tuning = TRUE, nfolds = 10, nlambda.inc = 20,
#                  nlambda.lat = 20, verbose = FALSE, parallel = TRUE, measure.inc = "c")
#    lambda.lat[k]<-coxem.c.k$selected.lambda.lat
#  }
#  table(lambda.inc)
#  table(lambda.lat)

## ----args4--------------------------------------------------------------------
args(cv_curegmifs)

## ----print--------------------------------------------------------------------
print(fitem)

## -----------------------------------------------------------------------------
summary(fitem)

## -----------------------------------------------------------------------------
summary(fit.cv)

## ----plot---------------------------------------------------------------------
plot(fitem)

## ----plotAIC------------------------------------------------------------------
plot(fitem, type = "cAIC")

## -----------------------------------------------------------------------------
plot(fit.cv)

## -----------------------------------------------------------------------------
coef.cAIC <- coef(fitem, model.select = "cAIC")

## -----------------------------------------------------------------------------
coef.5 <- coef(fitem, model.select = 5)

## -----------------------------------------------------------------------------
names(coef.cAIC)
all.equal(coef.cAIC$rate, coef.5$rate)
all.equal(coef.cAIC$alpha, coef.5$alpha)
all.equal(coef.cAIC$b0, coef.5$b0)
all.equal(coef.cAIC$beta_inc, coef.5$beta_inc)
all.equal(coef.cAIC$beta_lat, coef.5$beta_lat)

## ----pred---------------------------------------------------------------------
train.predict <- predict(fitem, model.select = "cAIC")

## -----------------------------------------------------------------------------
p_group <- ifelse(train.predict$p.uncured < 0.50, "Cured", "Susceptible")

## -----------------------------------------------------------------------------
km.cured <- survfit(Surv(cryr, relapse.death) ~ p_group, data = amltrain) 

## ----echo = FALSE, fig = TRUE-------------------------------------------------
plot(km.cured, mark.time = TRUE, lty = c(1,2), xlab = "Time (years)", ylab = "Relapse-free survival") 
legend(c(.9, .1), legend = c("Cured", "Susceptible"), lty = c(1, 2), bty = "n")

## -----------------------------------------------------------------------------
km.suscept <- survfit(Surv(cryr, relapse.death) ~ train.predict$latency.risk, data = amltrain, subset = (p_group == "Susceptible")) 

## ----echo = FALSE, fig = TRUE-------------------------------------------------
plot(km.suscept, mark.time = TRUE, lty = c(1,2), xlab = "Time (years)", ylab = "Relapse-free survival") 
legend(c(.9, .1), legend = c("Higher risk", "Lower risk"), lty = c(1,2), bty = "n")

## ----testpred-----------------------------------------------------------------
test.predict <- predict(fitem, newdata = amltest, model.select = "cAIC")

## -----------------------------------------------------------------------------
test_p_group <- ifelse(test.predict$p.uncured < 0.50, "Cured", "Susceptible")

## -----------------------------------------------------------------------------
km.cured.test <- survfit(Surv(cryr, relapse.death) ~ test_p_group, data = amltest) 

## ----echo = FALSE, fig = TRUE-------------------------------------------------
plot(km.cured.test, mark.time = TRUE, lty = c(1, 2), xlab = "Time (years)", ylab = "Relapse-free survival") 
legend(c(.4, .1), legend = c("Cured", "Susceptible"), lty = c(1,2), bty = "n")

## -----------------------------------------------------------------------------
km.suscept.test <- survfit(Surv(cryr, relapse.death) ~ test.predict$latency.risk, data = amltest, subset = (test_p_group == "Susceptible")) 

## ----echo = FALSE, fig = TRUE-------------------------------------------------
plot(km.suscept.test, mark.time = TRUE, lty = c(1,2), xlab = "Time (years)", ylab = "Relapse-free survival") 
legend(c(.4, .1), legend = c("Higher risk", "Lower risk"), lty = c(1, 2), bty = "n")

## -----------------------------------------------------------------------------
AUC(fitem, model.select = "cAIC")
AUC(fitem, newdata = amltest, model.select = "cAIC")
concordance_mcm(fitem, model.select = "cAIC")
concordance_mcm(fitem, newdata = amltest, model.select = "cAIC")

