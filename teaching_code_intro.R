# Load necessary packages
library(margins)   # For outcome regression
library(MatchIt)   # For causal inference matching
library(WeightIt)  # For causal inference weighting (IPW)
library(AIPW)      # For causal inference DR
library(survey)    # For survey-weighted estimation
library(tidyverse)     # For data manipulation
library(SuperLearner)
# Load the Lalonde dataset
data("lalonde", package = "Matching")

# View dataset structure
str(lalonde)

# ------------------- 1. Outcome Regression -------------------
options(scipen = 999)

# Fit Outcome Regression Model
outcome_model <- lm(re78 ~ treat + age + educ + black + hisp + married + nodegr + re74 + re75, data = lalonde)

summary(outcome_model)

# Compute ATE and 95% CI
marginal_effects <- margins(outcome_model, variables = "treat")
summary(marginal_effects)

# ------------------- 2. Propensity Score Estimation -------------------
ps_model <- glm(treat ~ age + educ + black + hisp + married + nodegr + re74 + re75, 
                family = binomial, data = lalonde)
lalonde$ps <- predict(ps_model, type = "response")

ggplot(lalonde, aes(x = ps, fill = factor(treat))) +
  geom_histogram(bins = 10, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("blue", "red"), labels = c("Control", "Treated")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +  # Set x-axis breaks at 0.1 intervals
  labs(title = "Propensity Score Distribution", x = "Propensity Score", y = "Count", fill = "Group") +
  theme_minimal()

# ------------------- 3. Propensity Score Matching (PSM) -------------------
match_out <- matchit(treat ~ age + educ + black + hisp + married + nodegr + re74 + re75, 
                     data = lalonde, method = "nearest")
matched_data <- match.data(match_out)

# Estimate treatment effect on matched data
matched_model <- lm(re78 ~ treat, data = matched_data)
summary(matched_model)

# ------------------- 4. Inverse Probability Weighting (IPW) using WeightIt -------------------
# Estimate weights using WeightIt package
ipw_weights <- weightit(treat ~ age + educ + black + hisp + married + nodegr + re74 + re75, 
                        data = lalonde, method = "ps")

# Apply IPW for treatment effect estimation
ipw_design <- svydesign(ids = ~1, weights = ~ipw_weights$weights, data = lalonde)
ipw_model <- svyglm(re78 ~ treat, design = ipw_design)
summary(ipw_model)

# ------------------- 5. Subclassification Based on Propensity Score -------------------
s.out1 <- matchit(treat ~ age + educ + black + hisp + married + nodegr + re74 + re75, data = lalonde,
                  method = "subclass", subclass = 7)
lalonde$subclass = s.out1$subclass
subclass_model = lm(re78 ~ treat*as.factor(subclass), data = lalonde)
marginal_effects <- margins(subclass_model, variables = "treat")
summary(marginal_effects)
# ------------------- 6. Doubly Robust Estimation -------------------
# Combines outcome regression and IPW
# Install and load AIPW package
# Estimate ATE using AIPW
cov = c("age", "educ", "black", "hisp", "married", "nodegr" , "re74" ,"re75")
W <- lalonde[, cov, drop=FALSE]
aipw_model <- AIPW$new(Y = lalonde$re78,
                       A = lalonde$treat,
                       W = W,
                       Q.SL.library = c("SL.lm"),
                       g.SL.library = c("SL.glm"),
                       k_split = 10,
                       verbose = FALSE)

# Extract results
aipw_model$fit()
aipw_model$summary()
print(aipw_model$result, digits = 2)


