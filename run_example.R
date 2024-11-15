library(netmeta)
library(nmajags)
library(tidyverse)
library(rlist)
library(rjags)
library(R2jags)


setwd("C:\\Users\\evrenogl\\Desktop\\Multiple_outcomes\\Data_and_Codes")

source("mvdata.R")

source("helpers.R")

source("jags_bivariate_NMA.R")

data <- read.csv("Data\\all_data_long.csv")

#data <- read.csv("Data\\dat2arm.csv")


dat <- mvdata(data = data,
              event = list(e.R,e.A),
              n = n,
              sm = "OR",
              studlab = study,
              treat = treat
)

reference.group <- "Placebo"

ref <- unname(which(dat$labtreat==reference.group))

run.data <- list(
  
  y = ifelse(is.na(dat$y),0,dat$y),
  var1 = dat$var$var1,
  var2 = dat$var$var2,
  ref = ref,
  Ns = dat$Ns,
  N2h = dat$N2h,
  T1 = dat$T[,1],
  T2 = dat$T[,2],
  T3 = dat$T[,3],
  NT = length(unique(dat$labtreat))
  
)

run = jags(
  data = run.data,
  inits = NULL,
  parameters.to.save = c(
    "res.ref1",
    "res.ref2",
    "rho1",
    "res.1",
    "res.2",
    "psi1",
    "psi2",
    "d1",
    "d2"
  ),
  n.chains = 2,
  n.iter = 10000,
  n.burnin = 2000,
  DIC = F,
  model.file = modfile
)


results = as.data.frame(run$BUGSoutput$summary)

results$ind = as.character(rownames(results))

####### Manipulate the results and create suitable datasets

#### Basic comparisons CA
basic_comp_1 = results %>%
  filter(grepl("ref1", ind))

row.names(basic_comp_1) <- dat$labtreat[-ref]

dat_treat1 <- dat$treat_out[[1]][-which(dat$treat_out[[1]]==reference.group)]

basic_comp_1 <- basic_comp_1[which(row.names(basic_comp_1) %in% dat_treat1),]

basic_comp_2 = results %>%
  filter(grepl("ref2", ind))

row.names(basic_comp_2) <- dat$labtreat[-ref]

dat_treat2 <- dat$treat_out[[2]][-which(dat$treat_out[[2]]==reference.group)]

basic_comp_2 <- basic_comp_2[which(row.names(basic_comp_2) %in% dat_treat2),]


## all results 
all_res1 <- results %>% 
  filter(grepl("res.1",ind))

## psi's (similar to tau)

psi <- results %>%
  filter(grepl(c("psi"), ind))

## rho

rho <- results %>% 
  filter(grepl("rho",ind))

### sims matrix 
sims_bugs_out <- as.data.frame(run$BUGSoutput$sims.matrix)

#### d's
d1 <- sims_bugs_out[,c(grepl("d1", names(sims_bugs_out)))]

d2 <- sims_bugs_out[,c(grepl("d2", names(sims_bugs_out)))]

names(d1) <- names(d2) <-  dat$labtreat

d1 <- d1 %>% 
  dplyr::select(all_of(dat$treat_out[[1]]))

d2 <- d2 %>% 
  dplyr::select(all_of(dat$treat_out[[2]]))

d1 <- as.matrix(d1)

d2 <- as.matrix(d2)

#rank1 <- rankogram(d1, pooled = "random", small.values = "undesirable")

#rank2 <- rankogram(d2, pooled = "random", small.values = "desirable")
       
# plot(rank1)
