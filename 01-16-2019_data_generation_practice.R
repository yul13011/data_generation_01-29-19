# This is the second practice coding for dissertation simulation study #
# In this coding, only one school level and one student level variable is simualted #
# Note all conditions are estimators are programmed #

# Set random seed #
# Step 1. Generate data #################################################################

# Step 1.1. Generate school level and student level covariates #
# number of schools = 1, ..., H, number of student per school = 1, ... ,K-h
# Set the number of schools to be 2000 #
# Generate the number of student per school to be N(200,80), min = 10 and max = 700
# Students per school is based on FL data
# Generate school-level V1 is continuous N(0, 1)
# Generate student level covariates x1 is continuous N(V1_h, 1)

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
seed <- 1234566
iter <- 1
if (identical(args, character(0))) {
    seed <- seed + iter
    message("No seed provided; using default seed ", seed)
} else {
    iter <- as.integer(args[[1]])
    seed <- seed + iter
    message("Setting seed to ", seed, " as requested")
}
set.seed(seed)

## Create invariant paramters.
params <- expand.grid(pi30 = 1,
                      pi31 = c(0.2, 2),
                      pi40 = c(0.2, 2),
                      pi41 = c(0.2, 2),
                      alpha0 = -2,
                      alpha1 = c(0.2, 2),
                      tau00 = -1,
                      tau01 = c(0.2, 2),
                      tau10 = c(0.2, 2),
                      tau11 = c(0.2, 2)) %>%
    ## Make it easier to interactively use.
    as_tibble() %>%
    ## Each row is a condition.
    mutate(condition = row_number()) %>%
    ## Move the condition column to be first.  "everything()" means
    ## all other columns.
    select(condition, everything())
## BEGIN FIXME
## For now subset the data so that we don't fill our RAM.
## Remove the stuff between the FIXME blocks when running this on the cluster.
params <- head(params, 3)
H <- 20 # 2000
## END FIXME
school <-
    tibble(schoolid = 1:H) %>%
    ## To merge schools with conditions, we need to repeat the schools
    ## as many times as we have conditions.
    replicate(nrow(params), ., simplify = FALSE) %>%
    bind_rows(.id = "condition") %>%
    ## Convert the condition column from character to integer for join.
    mutate(condition = as.integer(condition)) %>%
    inner_join(params, by = "condition") %>%
    ## All equations should apply across all schools.
    group_by(condition) %>%
    mutate(v1 = rnorm(n(), 0, 1), # v1~N(0,1)
           v2 = as.logical(rbinom(n(), 1, 0.5))) %>%
    mutate(phi0 = pi30 + pi31 * v1, # level 2 intercept as outcome of v1
           phi1 = pi40 + pi41 * v1, # level 2 slope as outcome of v1
           p2 = 1 / (1 + exp(-(alpha0 + alpha1 * v1))),
           s2 = as.logical(rbinom(n(), 1, p2)),
           z2 = as.logical(rbinom(n(), 1, 0.5)), # z2=1 indicates a sampled school was assigned the treatment condition
           eta0 = tau00 + tau01 * v1 + v2, # level 2 intercept as outcome of v1
           eta1 = tau10 + tau11 * v1 + v2) # level 2 slope as outcome of v1
## Predicted school selection probability.
fits <-
    school %>%
    ## Split data frames to be able to run glm.
    split(pull(school, condition)) %>%
    lapply(function(x) glm(s2 ~ v1, data = x,
                           family = binomial()))
prob_schl <- lapply(fits, predict,
                    type = "response") %>%
    unlist()
school <-
    school %>%
    bind_cols(prob_schl = prob_schl)
nstudents <- as.integer(rnorm(H,200,80))
nstudents[nstudents < 10] <- 10L
nstudents[nstudents > 700] <- 700L
data <-
    tibble(schoolid = rep.int(1:H, nstudents)) %>%
    inner_join(school, by = "schoolid") %>%
    ## All equations apply to all students in each school per condition.
    group_by(schoolid, condition) %>%
    ## Generate studentid from the group row number.
    mutate(studentid = row_number()) %>%
    ## n() is the group size.
    mutate(K = as.integer(rnorm(n(), 200, 80)),
           ## Truncate K to [10, 700]
           K = case_when(
               K < 10 ~ 10,
               K > 700 ~ 700)
           ) %>%
    mutate(x1 = rnorm(n(), v1, 1),
           y0 = v1 + x1 * v1, #potential outcome for treatment
           y1 = y0 + phi0 + phi1 * x1, # potential outcome for control
           TE = phi0 + phi1 * x1, # true individual treatment effect
           p1 = 1 / (1 + exp(-(eta0 + eta1 * x1))), # model for student selection
           s1 = as.logical(rbinom(n(), 1, p1)), # student selection follows a bernoulli distribution
           y = NA,
           y = case_when(s2 == 1 & s1 == 1 & z2 == 1 ~ y1,
                         s2 == 1 & s1 == 1 & z2 == 0 ~ y0))

file.name <- paste0("dataset", iter, ".dta")
foreign::write.dta(data, file.name)
