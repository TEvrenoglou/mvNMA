library(netmeta)
library(tidyverse)

setwd("C:\\Users\\evrenogl\\Desktop\\Multiple_outcomes\\Data_and_Codes")

data <- read.csv("Data\\all_data_long.csv")

#data <- read.csv("Data\\dat2arm.csv")

p1 <- pairwise(data = data,
              treat=treat,
              event = e.R,
              n = n,
              sm = "OR",
              studlab = study,
              append = F
              )

mod1 <- netmeta(p1,reference.group = "Placebo",common = F)

summary(mod1)

mod1$TE.random[,10]


p2 <- pairwise(data = data,
               treat=treat,
               event = e.A,
               n = n,
               sm = "OR",
               studlab = study
)


mod2 <- netmeta(p2,reference.group = "Placebo",common = F)

summary(mod2)

mod2$TE.random[,10]
