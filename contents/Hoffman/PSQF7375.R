library(readr)

library(lme4)
library(nlme)
library(car)

library(lmerTest) # type III anova, df adjustment (Kenward-Roger,...)

# Chapter/Example 6 -------------------------------------------------------
library()
ex6 <- read_csv("CHAPTER6.csv")
ex6$PersonID <- factor(ex6$PersonID)
ex6$session.f <- factor(ex6$session, levels = c(6, 1:5))

# Ch 6: 0: Saturated Means, Unstructured Variance Model
# TOTAL ANSWER KEY
summary(fit.gls <- gls(rt ~ session.f, data = ex6,
            correlation = corSymm(form = ~ 1 | PersonID)))
getVarCov(fit.gls)
cov2cor(getVarCov(fit.gls))

# Ch 6: 1b: Empty Means, Random Intercept Model
summary(fit <- lmer(rt ~ 1 + (1 | PersonID), data = ex6))
summary(fit.lme <- lme(rt ~ 1, random = ~ 1 | PersonID, data = ex6))
VarCorr(fit.lme)

fit.null <- gls(rt ~ 1, data = ex6) # null model
anova(fit.lme, fit.null)

# Ch 6: 2a: Fixed Linear Time, Random Intercept Model
summary(fit <- lmer(rt ~ time1 + (1 | PersonID), data = ex6))
summary(fit.lme <- lme(rt ~ time1, random = ~ 1 | PersonID, data = ex6))

anova(fit) # differs when using "lmerTest" package
anova.lme(fit.lme)
Anova(fit, test = "F") # refer to R companion

# Ch 6: 2b: Random Linear Time Model
summary(fit <- lmer(rt ~ time1 + (time1 | PersonID), data = ex6))
summary(fit.lme <- lme(rt ~ time1, random = ~ time1 | PersonID, data = ex6))

anova(fit)

library(emmeans)
emmeans::emmeans(fit, ~time1, at=list(time1 = c(0:5)))

# Ch 6: 3a: Fixed Quadratic, Random Linear Time Model
summary(fit <- lmer(rt ~ time1 + I(time1*time1) +
                            (time1 | PersonID), data = ex6))
summary(fit.lme <- lme(rt ~ time1 + I(time1*time1),
            random = ~ time1 | PersonID, data = ex6))
anova(fit)

fit.null <- gls(rt ~ time1 + I(time1*time1), data = ex6)
anova(fit.lme, fit.null) # null model has the same fixed effects, but no random effects at all

emmeans(fit, ~time1, at=list(time1 = c(0:5)))
emtrends(fit.lme, ~ time1, var="time1", at=list(time1 = c(0:5)))

# Ch 6: 3b: Random Quadratic Time Model
summary(fit <- lmer(rt ~ time1 + I(time1*time1) +
                            (time1 + I(time1*time1) | PersonID), data = ex6))
summary(fit.lme <- lme(rt ~ time1 + I(time1*time1),
                       random = ~ time1 + I(time1*time1) | PersonID, data = ex6))
Anova(fit, test = "F")

# Bonus: Testing AR1 Residual Correlation Across Time
summary(fit.lme <- lme(rt ~ time1 + I(time1*time1),
                       random = ~ time1 + I(time1*time1) 
                       | PersonID, 
                       data = ex6,
                       correlation = corAR1(form = ~ 1 | PersonID)))

mSt <- fit.lme$modelStruct
cSt <- mSt$corStruct
corMatrix(cSt)
summary(cSt)

cl <- getCall(fit.lme)
cl$correlation
intervals(fit.lme, which = "var-cov")$corStruct
getVarCov(fit.lme)
resid(fit.lme, type = "normalized")

# -------------------------------------------------------------------------


# Example 7 ---------------------------------------------------------------
ex7 <- read_csv("EXAMPLE7.csv")
ex7$ID <- factor(ex7$ID)

# TITLE1 "1a: Piecewise Unconditional Model - Random Early/Later Practice Slopes"

mod1 <- lmer(nm3rt ~ Slope12 + Slope26 + 
                     (Slope12 + Slope26 | ID), data = ex7)

mod1 <- lmer(nm3rt ~ time1 + I(time1^2) +
                     (time1 | ID) + (I(time1^2) | ID), data = ex7)
S(mod1)

mod2 <- lme(nm3rt ~ time1 + I(time1^2),
            random = ~ time1 + I(time1^2) | ID, data = ex7)
S(mod2)
