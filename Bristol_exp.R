#.libPaths()

suppressPackageStartupMessages({
  library(glmmTMB)
  library(car)
  library(popbio)
  library(ggplot2)
  library(ggbeeswarm)
  library(shades)
  library(RColorBrewer)
  library(DHARMa)
  library(Rmisc)
  library(dabestr)
  library(magrittr)
  library(tidyr)
  library(dplyr)
  library(ggeffects)
  library(lme4)
  library(optimx)
  library(nloptr)
  library(dfoptim)
  library(data.table)
  library(DHARMa)
  library(emmeans)
  library(ggpubr)
  library(patchwork)
  library(lubridate)
  library(performance)
  library(dplyr)#to check overdispersion in glmmTMB models
  library(survival)
  library(coxme) #for cox model
  library(purrr) #for forestplot
  library(ggpubr) #for forestplot
  library(survminer) #for forestplot
  library(AICcmodavg)
  library(forcats)
  library(fields) #for thin plate spline plots
  library(usethis)
  library(gitcreds)
  library(httr2)
  library(readr)
})



#edit_git_config() #put name and email in
#use_git()#this adds git repo
#create_github_token() #takes you to github to create personal access token
#gitcreds_set() #grant access to GitHub for RStudio
use_github()

R.version

## Experiment 1: Lifespan and reproduction for individuals in 10-15C temperatures either during development (egg stage to late L4, 'dev') or adulthood (24 hours from late L4 stage, 'adult')

#Reproduction data
adult <- read.csv('Bristol_adult.csv')
dev <- read.csv('Bristol_dev.csv')

#Lifespan data
ls <- read.csv('Bristol_lifespan.csv')


#Change from wide to long
adult_long <- adult %>% 
  pivot_longer(
    cols = `D1`:`D10`, 
    names_to = "Day",
    values_to = "value"
  )

adult_long <- na.omit(adult_long)

#Change from wide to long
dev_long <- dev %>% 
  pivot_longer(
    cols = `D1`:`D10`, 
    names_to = "Day",
    values_to = "value"
  )



