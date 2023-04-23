library(readxl)
library(lme4)

library(lmerTest)
library(tidyverse)
filename <- read_excel("BWPilot_CombinedData_20210803_fordryad_addvars_cleaned_noround2_v3.xlsx")

filename[filename=='.'] <- NA

data_set <- filename[, c(1, 2, 3, 4, 6)]

data_set <- data_set %>% drop_na()

data_set$`Post PANAS Positive Affect` <- as.numeric(data_set$`Post PANAS Positive Affect`)
data_set$`Pre PANAS Positive Affect` <- as.numeric(data_set$`Pre PANAS Positive Affect`)
data_set$`Round 1 Exercise` <- as.factor(data_set$`Round 1 Exercise`)
data_set$`Days from  Round1 Day1` <- as.numeric(data_set$`Days from  Round1 Day1`)


data_set$change <- data_set$`Post PANAS Positive Affect` - data_set$`Pre PANAS Positive Affect`

model = lmer(change ~ `Days from  Round1 Day1` + `Round 1 Exercise` * `Days from  Round1 Day1` +
               `Round 1 Exercise` + (1|ScrSubjectID),
             data=data_set,
             REML=TRUE)
summary(model)

anova(model)


