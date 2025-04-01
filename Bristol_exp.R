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
  library(Hmisc)
})



#edit_git_config() #put name and email in
#use_git()#this adds git repo
#create_github_token() #takes you to github to create personal access token
#gitcreds_set() #grant access to GitHub for RStudio
#use_github()

R.version

## Experiment 1: Lifespan and reproduction for individuals in 10-15C temperatures either during development (egg stage to late L4, 'dev') or adulthood (24 hours from late L4 stage, 'adult')

#Reproduction data
adult <- read.csv('Bristol_adult.csv')
dev <- read.csv('Bristol_dev.csv')

#Lifespan data
ls <- read.csv('Bristol_lifespan.csv')

#Split lifespan data into developmental and adulthood daf-2 treatments
dev_ls <- subset(ls, Treatment == 'dev')
adult_ls <- subset(ls, Treatment == 'adult')

#Lifespan data including mutants
mut_LS <- read.csv('Bristol_LS.csv')


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

adult_long$value <- as.numeric(adult_long$value)
adult_long$Day <- as.factor(adult_long$Day)

dev_long$value <- as.numeric(dev_long$value)
dev_long$Day <- as.factor(dev_long$Day)


ad_long_daf <- subset(adult_long, Treatment == 'daf')
ad_long_ev <- subset(adult_long, Treatment == 'ev')


level_order <- c('D1', 'D2', 'D3', 'D4','D5','D6','D7','D8','D9','D10')



##### Plots for adult daf-2 reproduction

## Age-specific reproduction plot

adult_rep <-ggplot(data=adult_long, aes(x=factor(Day, levels = level_order), y=value, group=Treatment, color=Treatment))+
  geom_jitter(alpha = 0.2, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5))+
  stat_summary(fun.data="mean_cl_boot", geom="errorbar", size = 1, width=0.0, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="point", size = 3, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="line",  size=1, position = position_dodge(0.5)) +
  theme_classic()+
  labs(y="Offspring number", x="", size=20)+
  labs(col="")+
  theme(axis.title.y = element_text(size=12))+
  theme(axis.title.x = element_text(size=12))+
  scale_color_manual(values=c('deepskyblue4','darkorange2'))+
  theme(legend.key.width = unit(0.5,"cm"))+
  coord_cartesian(ylim = c(0,85))+
  theme(legend.position = c(0.8,0.9))


#Early reproduction for daf-2 worms looks higher
adult_rep

ggsave('adult_rep_plot.pdf', height = 20, width = 25, units = 'cm')

## Calculate lambda

L <- matrix(nrow = nrow(adult), ncol = 3)


for (i in 1:nrow(adult)){
  Les <- matrix(0, ncol = 12, nrow = 12)
  diag(Les[-1,]) <- rep(1, 11) # add the 1s for survival probability diagonally
  Fert <- c(0,0, as.numeric(as.vector(adult[i,][4:13]))) #the columns in data that has the reproductive counts
  Fert[is.na(Fert)] <- 0 #makes all NAs into 0s
  Les[1,] <- c(Fert)
  class(Les) <- "leslie.matrix"
  Lambda <- popbio::eigen.analysis(Les)$lambda1
  L[i, 1:3] <- c(paste0(adult$Treatment[i]), paste0(adult$ID[i]), Lambda)
  
}

colnames(L)<-c("Treatment", "ID", "Lambda")

Data<-as.data.frame(L)

Data$Lambda<-as.numeric(as.character(Data$Lambda))

Lambda_adult <-
  Data %>%
  load(Treatment, Lambda,
         idx = list(c("ev", "daf")))

vector<-c('deepskyblue4','darkorange2')

## Lambda dabestr plot
Lambda_plot_adult <- dabest_plot(mean_diff(Lambda_adult),
                          palette = c('deepskyblue4','darkorange2'), axes.title.fontsize=13, 
                          rawplot.ylabel="Fitness", effsize.ylabel = "", rawplot.markersize=2, tick.fontsize=9)
#Daf-2 worms have higher fitness in natural temperatures
Lambda_plot_adult

ggsave('lambda_adult.pdf', height = 20, width = 25, units = 'cm' )


## LRS
totalrep<-na.omit(as.data.frame.table(tapply(adult_long$value,list(adult_long$Treatment, adult_long$ID),sum)))

names(totalrep)<-c("Treatment", "Replicate", "Totrep")

Totrep <-
  totalrep %>%
  load(Treatment, Totrep,
         idx = c('ev','daf'))

LRS <- dabest_plot(mean_diff(Totrep),
            #color.column = Treatment, #this remove the legend - make sure colours are correct though
            palette = c('darkorange2','deepskyblue4'), axes.title.fontsize=13, 
            rawplot.ylabel="Total reproduction", effsize.ylabel = "Mean difference", rawplot.markersize=2, tick.fontsize=9)

LRS

adult_repro <- ggarrange(adult_rep, 
                         ggarrange(Lambda_plot_adult, LRS, ncol = 2),
                         nrow = 2,
                         common.legend = TRUE)
adult_repro

ggsave('all_adult_rep.pdf', plot = adult_repro,
       width = 10, height = 10, units = 'in')

