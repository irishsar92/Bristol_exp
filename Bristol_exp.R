


# Setup 
```{r}
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
adult <- read.csv('Bristol_adult_rep.csv')
dev <- read.csv('Bristol_dev_rep.csv')

#Lifespan data
ls <- read.csv('Bristol_lifespan.csv')

#Split lifespan data into developmental and adulthood daf-2 treatments
dev_ls <- subset(ls, Treatment == 'dev')
adult_ls <- subset(ls, Treatment == 'adult')

#Lifespan data including mutants
mut_LS <- read.csv('Bristol_LS_mut.csv')
mut_rep <- read.csv('Bristol_repro_mut.csv')


level_order <- c('D1', 'D2', 'D3', 'D4','D5','D6','D7','D8','D9','D10')

```

# N2 only adult and dev daf-2 RNAi in 10 to 15C daily gradual thermocycle - format data
```{r}
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

```



##### Plots for adult daf-2 reproduction

```{r}
## Age-specific reproduction plot

adult_rep <-ggplot(data=adult_long, aes(x=factor(Day, levels = level_order), y=value, group=Treatment, color=Treatment))+
  geom_jitter(alpha = 0.2, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5))+
  stat_summary(fun.data="mean_cl_boot", geom="errorbar", size = 1, width=0.0, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="point", size = 3, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="line",  size=1, position = position_dodge(0.5)) +
  theme_classic()+
  labs(y="Offspring number", x="", size=20)+
  labs(col="")+
  theme(axis.title.y = element_text(size=16))+
  theme(axis.title.x = element_text(size=14))+
  scale_color_manual(values=c('#4DBBD5FF','#E64B35FF'))+
  theme(legend.key.width = unit(0.5,"cm"))+
  coord_cartesian(ylim = c(0,85))+
  theme(legend.position = c(0.8,0.9))


#Early reproduction for daf-2 worms looks higher
adult_rep

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

Lambda_dab <- mean_diff(Lambda_adult)

## Lambda dabestr plot
Lambda_plot_adult <- dabest_plot(Lambda_dab, FALSE, swarm_label = 'Lambda', raw_marker_spread = 1, custom_palette = 'npg', swarm_x_text = 12, swarm_y_text = 16, contrast_y_text = 16, contrast_x_text = 12, raw_marker_alpha = 0.3, tufte_size = 1)
Lambda_plot_adult


## LRS
totalrep<-na.omit(as.data.frame.table(tapply(adult_long$value,list(adult_long$Treatment, adult_long$ID),sum)))

names(totalrep)<-c("Treatment", "Replicate", "Totrep")

Totrep <-
  totalrep %>%
  load(Treatment, Totrep,
         idx = c('ev','daf'))

Totrep_dab <- mean_diff(Totrep)

LRS_adult <- dabest_plot(Totrep_dab, FALSE, swarm_label = 'LRS', raw_marker_spread = 1, custom_palette = 'npg', swarm_x_text = 12, swarm_y_text = 16, contrast_y_text = 16, contrast_x_text = 12, raw_marker_alpha = 0.3, tufte_size = 1)

LRS_adult

adult_repro <- ggarrange(adult_rep, 
                         ggarrange(Lambda_plot_adult, LRS_adult, ncol = 2),
                         nrow = 2,
                         common.legend = TRUE)
adult_repro

ggsave('all_adult_rep.pdf', plot = adult_repro,
       width = 10, height = 10, units = 'in')

```


# Developmental daf-2 RNAi reproduction

```{r}
## Age-specific reproduction plots

 dev_rep <-ggplot(data=dev_long, aes(x=factor(Day, levels = level_order), y=value, group=Treatment, color=Treatment))+
  geom_jitter(alpha = 0.2, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5))+
  stat_summary(fun.data="mean_cl_boot", geom="errorbar", size = 1, width=0.0, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="point", size = 3, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="line",  size=1, position = position_dodge(0.5)) +
  theme_classic()+
  labs(y="Offspring number", x="Day", size=20)+
  labs(col="")+
  theme(axis.title.y = element_text(size=16))+
  theme(axis.title.x = element_text(size=14))+
  scale_color_manual(values=c('#4DBBD5FF','#E64B35FF'))+
  theme(legend.key.width = unit(0.5,"cm"))+
  coord_cartesian(ylim = c(0,85))+
  theme(legend.position = c(0.8,0.9))

dev_rep


## Calculate lambda

L <- matrix(nrow = nrow(dev), ncol = 3)

for (i in 1:nrow(dev)){
  Les <- matrix(0, ncol = 12, nrow = 12)
  diag(Les[-1,]) <- rep(1, 11) # add the 1s for survival probability diagonally
  Fert <- c(0,0, as.numeric(as.vector(dev[i,][4:13]))) #the columns in data that has the reproductive counts
  Fert[is.na(Fert)] <- 0 #makes all NAs into 0s
  Les[1,] <- c(Fert)
  class(Les) <- "leslie.matrix"
  Lambda <- popbio::eigen.analysis(Les)$lambda1
  L[i, 1:3] <- c(paste0(dev$Treatment[i]), paste0(dev$ID[i]), Lambda)
  
}

print(Fert)
colnames(L)<-c("Treatment", "ID", "Lambda")

Data_dev<-as.data.frame(L)

Data_dev$Lambda<-as.numeric(as.character(Data_dev$Lambda))

Lambda_dev <-
  Data_dev %>%
  load(Treatment, Lambda,
         idx = list(c("ev", "daf")))

Lambda_dab <- mean_diff(Lambda_dev)

## Lambda dabestr plot
Lambda_plot_dev <- dabest_plot(Lambda_dab, FALSE, swarm_label = 'Lambda', raw_marker_spread = 1, custom_palette = 'npg', swarm_x_text = 12, swarm_y_text = 16, contrast_y_text = 16, contrast_x_text = 12, raw_marker_alpha = 0.3, tufte_size = 1)


Lambda_plot_dev


## LRS

totalrep<-na.omit(as.data.frame.table(tapply(dev_long$value,list(dev_long$Treatment, dev_long$ID),sum)))

names(totalrep)<-c("Treatment", "Replicate", "Totrep")

Totrep <-
  totalrep %>%
  load(Treatment, Totrep,
         idx = c('ev','daf'))

Totrep_dab <- mean_diff(Totrep)

LRS_dev <- dabest_plot(Totrep_dab, FALSE, swarm_label = 'LRS', raw_marker_spread = 1, custom_palette = 'npg', swarm_x_text = 12, swarm_y_text = 16, contrast_y_text = 16, contrast_x_text = 12, raw_marker_alpha = 0.3, tufte_size = 1)
LRS_dev

dev_repro <- ggarrange(dev_rep, 
          ggarrange(Lambda_plot_dev, LRS_dev, ncol = 2),
          nrow = 2,
          common.legend = TRUE)
dev_repro

ggsave('all_dev_rep.pdf', plot = dev_repro,
       width = 10, height = 10, units = 'in')

```


# Adult daf-2 lifespan
```{r}

class(adult_ls$Treatment_ID)
adult_ls$Treatment_ID <- as.factor(adult_ls$Treatment_ID)
adult_ls$Treatment_ID <- relevel(adult_ls$Treatment_ID, "ev")

adult_ls$Cause <- as.factor(adult_ls$Cause)

levels(adult_ls$Cause)

 palette <- c('#E64B35FF','#4DBBD5FF')

adult_surv <- survfit(Surv(Age,Event)~Treatment_ID,data=adult_ls)

adult_lifespan<-ggsurvplot(adult_surv, data = adult_ls, ylab="Survival probability", size= 0.8, 
                          font.ylab= 16, font.xlab= 16, palette=palette, legend = "none", legend.title = "", 
                          legend.labs = c("ev", "daf-2"), 
                          title = "", censor = FALSE, xlab = "Day", xlim=c(0,65), break.time.by = 10, 
                          position= position_dodge(0.5), font.tickslab = c(9))
adult_lifespan


adult_LS_matcen <- adult_ls %>%
  filter(Cause != 'M')

adult_surv_matcen <- survfit(Surv(Age,Event)~Treatment_ID,data=adult_LS_matcen)

adult_lifespan_matcen<-ggsurvplot(adult_surv_matcen, data = adult_ls, ylab="Survival probability", size= 0.8, 
                          font.ylab= 16, font.xlab= 16, palette=palette, legend = "none", legend.title = "", 
                          legend.labs = c("ev", "daf-2"), 
                          title = "", censor = FALSE, xlab = "Day", xlim=c(0,65), break.time.by = 10,
                          position= position_dodge(0.5), font.tickslab = c(12))

adult_lifespan_matcen


meforest <- function(cox, ev){  #Eds version of forest plot
  require(AICcmodavg)
  require(ggplot2)
  require(forcats)
  store <- matrix(nrow = length(cox$coefficients) + 1, ncol = 4)
  ref <- c(paste(ev), 0, NA, NA)
  store[1,] <- ref
  for (x in 1:length(cox$coefficients)){
    y = x+1
    mean<-cox$coefficients[x]
    emean <- (mean)
    CIL <- cox$coefficients[x] - (1.96*extractSE(cox)[x])
    CUL <- cox$coefficients[x] + (1.96*extractSE(cox)[x])
    eCIL <- (CIL)
    eCUL <- (CUL)
    store[y, 1:4] <- c(names(cox$coefficients[x]), emean, eCIL, eCUL)}
  colnames(store) <- c("Treatment_ID", "mean", "CIL","CUL")
  store <- as.data.frame(store)
  store$mean <- as.numeric(as.character(store$mean))
  store$CIL <- as.numeric(as.character(store$CIL))
  store$CUL <- as.numeric(as.character(store$CUL))
  forest<-ggplot(store, aes(x = mean, xmax = CUL, xmin = CIL, y = Treatment_ID, colour = Treatment_ID)) +
    geom_point(size=2, shape = 19, colour = c('#E64B35FF','#4DBBD5FF')) +
    geom_errorbarh(height=0.25, size=0.8,  colour = c('#E64B35FF','#4DBBD5FF')) +
    theme_minimal() +
    geom_vline(xintercept=0, linetype="dotted", size=0.8) +
    xlab("Hazard Ratio")+
    ylab("")+
      theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 16))+
    scale_y_discrete(limits = (levels(store$Treatment_ID)), labels = c('ev','daf-2'))+
    expand_limits(x = c(-1.5,0.5))
  return(forest)
}

cox_LS <- coxme(Surv(Age, Event) ~ Treatment_ID + (1|Plate_ID), data = adult_ls)

ForestAll<-meforest(cox_LS, "ev")

ForestAll

cox_matcen <- coxme(Surv(Age, Event) ~ Treatment_ID + (1|Plate_ID), data = adult_LS_matcen)

Forest_matcen <- meforest(cox_matcen, 'ev')
Forest_matcen



LS_all_plot <- ggarrange(adult_lifespan$plot, adult_lifespan_matcen$plot, ForestAll, Forest_matcen, nrow = 2, ncol = 2,heights = c(2, 1), labels = c('A', 'B'))
LS_all_plot
ggsave('BristolLS_adult_plot.pdf', height = 8, width = 12)
```


# Dev daf-2 LS

```{r}

class(dev_ls$Treatment_ID)
dev_ls$Treatment_ID <- as.factor(dev_ls$Treatment_ID)
dev_ls$Treatment_ID <- relevel(dev_ls$Treatment_ID, "ev")

dev_ls$Cause <- as.factor(dev_ls$Cause)

levels(dev_ls$Cause)

 palette <- c('#E64B35FF','#4DBBD5FF')

dev_surv <- survfit(Surv(Age,Event)~Treatment_ID,data=dev_ls)

dev_lifespan<-ggsurvplot(dev_surv, data = dev_ls, ylab="Survival probability", size= 0.8, 
                          font.ylab= 16, font.xlab= 16, palette=palette, legend = "none", legend.title = "", 
                          legend.labs = c("ev", "daf-2"), 
                          title = "", censor = FALSE, xlab = "Day", xlim=c(0,65), break.time.by = 10, 
                          position= position_dodge(0.5), font.tickslab = c(9))
dev_lifespan


dev_LS_matcen <- dev_ls %>%
  filter(Cause != 'M')

dev_surv_matcen <- survfit(Surv(Age,Event)~Treatment_ID,data=dev_LS_matcen)

dev_lifespan_matcen<-ggsurvplot(dev_surv_matcen, data = dev_ls, ylab="Survival probability", size= 0.8, 
                          font.ylab= 16, font.xlab= 16, palette=palette, legend = "none", legend.title = "", 
                          legend.labs = c("ev", "daf-2"), 
                          title = "", censor = FALSE, xlab = "Day", xlim=c(0,65), break.time.by = 10,
                          position= position_dodge(0.5), font.tickslab = c(12))

dev_lifespan_matcen


meforest <- function(cox, ev){  #Eds version of forest plot
  require(AICcmodavg)
  require(ggplot2)
  require(forcats)
  store <- matrix(nrow = length(cox$coefficients) + 1, ncol = 4)
  ref <- c(paste(ev), 0, NA, NA)
  store[1,] <- ref
  for (x in 1:length(cox$coefficients)){
    y = x+1
    mean<-cox$coefficients[x]
    emean <- (mean)
    CIL <- cox$coefficients[x] - (1.96*extractSE(cox)[x])
    CUL <- cox$coefficients[x] + (1.96*extractSE(cox)[x])
    eCIL <- (CIL)
    eCUL <- (CUL)
    store[y, 1:4] <- c(names(cox$coefficients[x]), emean, eCIL, eCUL)}
  colnames(store) <- c("Treatment_ID", "mean", "CIL","CUL")
  store <- as.data.frame(store)
  store$mean <- as.numeric(as.character(store$mean))
  store$CIL <- as.numeric(as.character(store$CIL))
  store$CUL <- as.numeric(as.character(store$CUL))
  forest<-ggplot(store, aes(x = mean, xmax = CUL, xmin = CIL, y = Treatment_ID, colour = Treatment_ID)) +
    geom_point(size=2, shape = 19, colour = c('#E64B35FF','#4DBBD5FF')) +
    geom_errorbarh(height=0.25, size=0.8,  colour = c('#E64B35FF','#4DBBD5FF')) +
    theme_minimal() +
    geom_vline(xintercept=0, linetype="dotted", size=0.8) +
    xlab("Hazard Ratio")+
    ylab("")+
      theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 16))+
    scale_y_discrete(limits = (levels(store$Treatment_ID)), labels = c('ev','daf-2'))+
    expand_limits(x = c(-1.5,0.5))
  return(forest)
}

cox_LS <- coxme(Surv(Age, Event) ~ Treatment_ID + (1|Plate_ID), data = dev_ls)

ForestAll<-meforest(cox_LS, "ev")

ForestAll

cox_matcen <- coxme(Surv(Age, Event) ~ Treatment_ID + (1|Plate_ID), data = dev_LS_matcen)

Forest_matcen <- meforest(cox_matcen, 'ev')
Forest_matcen



LS_dev_all_plot <- ggarrange(dev_lifespan$plot, dev_lifespan_matcen$plot, ForestAll, Forest_matcen, nrow = 2, ncol = 2,heights = c(2, 1), labels = c('A', 'B'))
LS_all_plot
ggsave('BristolLS_dev_plot.pdf', height = 8, width = 12)

```


# N2s + SKN-1 and DAF-16 Mutants
```{r}

repro_mut <- read.csv('Bristol_repro_mut.csv')

mut_long <- repro_mut %>% 
  pivot_longer(
    cols = `D1`:`D15`, 
    names_to = "Day",
    values_to = "value"
)

mut_long$Day <- as.factor(mut_long$Day)
str(mut_long)

levels(mut_long$Day) <- list('1' = 'D1', '2' = 'D2', '3' = 'D3', '4' = 'D4', '5' = 'D5', '6' = 'D6', '7' = 'D7', '8' = 'D8', '9' = 'D9', '10' = 'D10', '11' = 'D11', '12' = 'D12', '13' = 'D13', '14' = 'D14', '15' = 'D15')

mut_long$Day <- as.factor(mut_long$Day)
#mut_long <- na.omit(mut_long)
str(mut_long)

SKN <- mut_long %>% 
  filter(Strain == 'SKN')

DAF16 <- mut_long %>% 
  filter(Strain == 'DAF16')

N2 <- mut_long %>%
  filter(Strain == 'N2')

mut_long <- mut_long %>%
  unite(Tr_str, c('Treatment','Strain'), remove= F)

mut_long$Tr_str <- factor(mut_long$Tr_str, levels = c("ev_N2", "daf_N2", "ev_SKN", "daf_SKN"))


totalrep<-na.omit(as.data.frame.table(tapply(mut_long$value,list(mut_long$Tr_str, mut_long$ID),sum)))



names(totalrep)<-c("Treatment", "Replicate", "Totrep")


Totrep <-
  totalrep %>%
  load(Treatment, Totrep,
         idx = c('ev_N2', 'daf_N2', 'ev_SKN','daf_SKN'))

LRS <- dabest_plot(mean_diff(Totrep),
            #color.column = Treatment, #this remove the legend - make sure colours are correct though
         axes.title.fontsize=13, 
            rawplot.ylabel="Total reproduction", effsize.ylabel = "Mean difference", rawplot.markersize=2, tick.fontsize=9)

LRS
Totrep <-
  totalrep %>%
  dabest(Treatment, Totrep,
         idx = c('ev_N2','daf_N2', 'ev_SKN', 'daf_SKN','ev_DAF16','daf_DAF16'))

LRS <- plot(mean_diff(Totrep),
          #color.column = Treatment, #this remove the legend - make sure colours are correct though
          axes.title.fontsize=13, 
          rawplot.ylabel="Total reproduction", effsize.ylabel = "Mean difference", rawplot.markersize=2, tick.fontsize=9)

LRS

 SKN_rep <-ggplot(data=SKN, aes(x=factor(Day), y=value, group=Treatment, color=Treatment))+
  geom_jitter(alpha = 0.2, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5))+
  stat_summary(fun.data="mean_cl_boot", geom="errorbar", size = 1, width=0.0, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="point", size = 3, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="line",  size=1, position = position_dodge(0.5)) +
  theme_classic()+
  labs(y="Offspring number", x="Day", size=20)+
  labs(col="")+
  theme(axis.title.y = element_text(size=14))+
  theme(axis.title.x = element_text(size=14))+
  scale_color_manual(values=c('deepskyblue4','darkorange2'))+
  theme(legend.key.width = unit(0.5,"cm"))+
  coord_cartesian(ylim = c(0,85))+
  theme(legend.position = c(0.8,0.9))

SKN_rep



SKN_wide <- spread(SKN, Day, value)

within(SKN_wide, rm(Strain))
## Calculate lambda

L <- matrix(nrow = nrow(SKN_wide), ncol = 3)


for (i in 1:nrow(SKN_wide)){
  Les <- matrix(0, ncol = 9, nrow = 9)
  diag(Les[-1,]) <- rep(1, 8) # add the 1s for survival probability diagonally
  Fert <- c(0,0, as.numeric(as.vector(SKN_wide[i,][4:10]))) #the columns in data that has the reproductive counts
  Fert[is.na(Fert)] <- 0 #makes all NAs into 0
  Les[1,] <- c(Fert)
  class(Les) <- "leslie.matrix"
  Lambda <- popbio::eigen.analysis(Les)$lambda1
  L[i, 1:3] <- c(paste0(SKN_wide$Treatment[i]), paste0(SKN_wide$ID[i]), Lambda)
  
}



colnames(L)<-c("Treatment", "ID", "Lambda")

Data<-as.data.frame(L)

Data$Lambda<-as.numeric(as.character(Data$Lambda))

Lambda_SKN <-
  Data %>%
  dabest(Treatment, Lambda,
         idx = list(c("ev", "daf")))

vector<-c('deepskyblue4','darkorange2')


## Lambda dabestr plot

Lambda_plot_SKN <- plot(mean_diff(Lambda_SKN),
     #color.column = Treatment,
     palette = c('darkorange2','deepskyblue4'), axes.title.fontsize=13, 
                 rawplot.ylabel="Fitness", effsize.ylabel = "", rawplot.markersize=2, tick.fontsize=9)

Lambda_plot_SKN

 DAF16_rep <-ggplot(data=DAF16, aes(x=factor(Day), y=value, group=Treatment, color=Treatment))+
  geom_jitter(alpha = 0.2, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5))+
  stat_summary(fun.data="mean_cl_boot", geom="errorbar", size = 1, width=0.0, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="point", size = 3, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="line",  size=1, position = position_dodge(0.5)) +
  theme_classic()+
  labs(y="Offspring number", x="Day", size=20)+
  labs(col="")+
  theme(axis.title.y = element_text(size=14))+
  theme(axis.title.x = element_text(size=14))+
  scale_color_manual(values=c('deepskyblue4','darkorange2'))+
  theme(legend.key.width = unit(0.5,"cm"))+
  coord_cartesian(ylim = c(0,85))+
  theme(legend.position = c(0.8,0.9))
 
 
 
 
 

DAF16_rep


 N2_rep <-ggplot(data=N2, aes(x=factor(Day), y=value, group=Treatment, color=Treatment))+
  geom_jitter(alpha = 0.2, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5))+
  stat_summary(fun.data="mean_cl_boot", geom="errorbar", size = 1, width=0.0, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="point", size = 3, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="line",  size=1, position = position_dodge(0.5)) +
  theme_classic()+
  labs(y="Offspring number", x="Day", size=20)+
  labs(col="")+
  theme(axis.title.y = element_text(size=14))+
  theme(axis.title.x = element_text(size=14))+
  scale_color_manual(values=c('deepskyblue4','darkorange2'))+
  theme(legend.key.width = unit(0.5,"cm"))+
  coord_cartesian(ylim = c(0,85))+
  theme(legend.position = c(0.8,0.9))
 
 N2_rep
 

DAF16_wide <- spread(DAF16, Day, value)

within(DAF16_wide, rm(Strain))
## Calculate lambda

L <- matrix(nrow = nrow(DAF16_wide), ncol = 3)


for (i in 1:nrow(DAF16_wide)){
  Les <- matrix(0, ncol = 9, nrow = 9)
  diag(Les[-1,]) <- rep(1, 8) # add the 1s for survival probability diagonally
  Fert <- c(0,0, as.numeric(as.vector(DAF16_wide[i,][4:10]))) #the columns in data that has the reproductive counts
  Fert[is.na(Fert)] <- 0 #makes all NAs into 0
  Les[1,] <- c(Fert)
  class(Les) <- "leslie.matrix"
  Lambda <- popbio::eigen.analysis(Les)$lambda1
  L[i, 1:3] <- c(paste0(DAF16_wide$Treatment[i]), paste0(DAF16_wide$ID[i]), Lambda)
  
}



colnames(L)<-c("Treatment", "ID", "Lambda")

Data<-as.data.frame(L)

Data$Lambda<-as.numeric(as.character(Data$Lambda))

Lambda_DAF16 <-
  Data %>%
  dabest(Treatment, Lambda,
         idx = list(c("ev", "daf")))

vector<-c('deepskyblue4','darkorange2')


## Lambda dabestr plot

Lambda_plot_DAF16 <- plot(mean_diff(Lambda_DAF16),
     #color.column = Treatment,
     palette = c('darkorange2','deepskyblue4'), axes.title.fontsize=13, 
                 rawplot.ylabel="Fitness", effsize.ylabel = "", rawplot.markersize=2, tick.fontsize=9)

Lambda_plot_DAF16
 
 


rep <- ggplot(data=mut_long, aes(x=factor(Day), y=value, group = Tr_str, color=Strain))+
  geom_jitter(alpha = 0.2, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6))+
  stat_summary(fun.data="mean_cl_boot", geom="errorbar", size = 1, width=0.0, position = position_dodge(0.6)) +
  stat_summary(fun.data="mean_cl_boot", geom="point", size = 3, position = position_dodge(0.6)) +
  stat_summary(fun.data="mean_cl_boot", geom="line",  size=1, position = position_dodge(0.6), aes(linetype = Treatment)) +
  theme_classic()+
  labs(y="Offspring number", x="Day", size=20)+
  labs(col="")+
  theme(axis.title.y = element_text(size=14))+
  theme(axis.title.x = element_text(size=14))+
  theme(legend.key.width = unit(0.5,"cm"))+
  scale_color_manual(values = c('deepskyblue4','darkorange2','goldenrod'))+
  coord_cartesian(ylim = c(0,85))+
  theme(legend.position = c(0.8,0.9))
rep



N2_wide <- spread(N2, Day, value)

within(N2_wide, rm(Strain))
## Calculate lambda

L <- matrix(nrow = nrow(N2_wide), ncol = 3)


for (i in 1:nrow(N2_wide)){
  Les <- matrix(0, ncol = 9, nrow = 9)
  diag(Les[-1,]) <- rep(1, 8) # add the 1s for survival probability diagonally
  Fert <- c(0,0, as.numeric(as.vector(N2_wide[i,][4:10]))) #the columns in data that has the reproductive counts
  Fert[is.na(Fert)] <- 0 #makes all NAs into 0s
  Les[1,] <- c(Fert)
  class(Les) <- "leslie.matrix"
  Lambda <- popbio::eigen.analysis(Les)$lambda1
  L[i, 1:3] <- c(paste0(N2_wide$Treatment[i]), paste0(N2_wide$ID[i]), Lambda)
  
}



colnames(L)<-c("Treatment", "ID", "Lambda")

Data<-as.data.frame(L)

Data$Lambda<-as.numeric(as.character(Data$Lambda))

Lambda_N2 <-
  Data %>%
  dabest(Treatment, Lambda,
         idx = list(c("ev", "daf")))

Lambda_N2

vector<-c('deepskyblue4','darkorange2')


## Lambda dabestr plot

Lambda_plot_N2 <- plot(mean_diff(Lambda_N2),
     #color.column = Treatment,
     palette = c('darkorange2','deepskyblue4'), axes.title.fontsize=13, 
                 rawplot.ylabel="Fitness", effsize.ylabel = "", rawplot.markersize=2, tick.fontsize=9)

Lambda_plot_N2 + Lambda_plot_DAF16+Lambda_plot_SKN


```


##Grouping datasets

```{r}
fin <- read.csv('Finn_data.csv')
repro <- read.csv('Bristol_adult.csv')
repro2 <- read.csv('Bristol_repro_mut.csv')

str(fin)
fin$infection <- NULL
fin$matracide <- NULL
fin$strain_treat <- NULL
str(fin)

str(repro)

str(repro2)
repro2$D11 <- NULL
repro2$D12 <- NULL
repro2$D13 <- NULL
repro2$D14 <- NULL
repro2$D15 <- NULL
str(repro2)


all<- rbind(fin, repro, repro2)
#Add individual ID - since ID numbers overlap across datasets
all <- tibble::rowid_to_column(all, "ID_ind")

all <- all%>%
  unite(Tr_str, c('Treatment','Strain'), remove= F)

```


##Plots for combined data

```{r}

#Get data in long form
all_long <- all %>% 
  pivot_longer(
    cols = `D1`:`D10`, 
    names_to = "Day",
    values_to = "value"
)

#Change Day to a factor so I can change levels to numeric values
all_long$Day <- as.factor(all_long$Day)

#Give levels for Day a numeric value
levels(all_long$Day) <- list('1' = 'D1', '2' = 'D2', '3' = 'D3', '4' = 'D4', '5' = 'D5', '6' = 'D6', '7' = 'D7', '8' = 'D8', '9' = 'D9', '10' = 'D10')

#Change class to numeric
all_long$Day <- as.numeric(all_long$Day)

#Get only N2s 
N2_all <- all_long %>%
  filter(Strain == 'N2')

###Need to exclude Lucas' data - not sure if he swapped mutants around
N2_all <- N2_all %>%
  filter(N2_all$Counter != 'lucas')

 N2_rep <-ggplot(data=N2_all, aes(x=factor(Day), y=value, group=Treatment, color=Treatment))+
  geom_jitter(alpha = 0.2, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5))+
  stat_summary(fun.data="mean_cl_boot", geom="errorbar", size = 1, width=0.0, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="point", size = 3, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="line",  size=1, position = position_dodge(0.5)) +
  theme_classic()+
  labs(y="Offspring number", x="Day", size=20)+
  labs(col="")+
  theme(axis.title.y = element_text(size=14))+
  theme(axis.title.x = element_text(size=14))+
  scale_color_manual(values=c('deepskyblue4','darkorange2'))+
  theme(legend.key.width = unit(0.5,"cm"))+
  coord_cartesian(ylim = c(0,85))+
  theme(legend.position = c(0.8,0.9))
 
 N2_rep
 
#Spread data back out
N2_wide <- spread(N2_all, Day, value)
within(N2_wide, rm(Strain))

N2_wide <- na.omit(N2_wide)

## Calculate lambda
L <- matrix(nrow = nrow(N2_wide), ncol = 3)

for (i in 1:nrow(N2_wide)){
  Les <- matrix(0, ncol = 9, nrow = 9)
  diag(Les[-1,]) <- rep(1, 8) # add the 1s for survival probability diagonally
  Fert <- c(0,0, as.numeric(as.vector(N2_wide[i,][4:10]))) #the columns in data that has the reproductive counts
  Fert[is.na(Fert)] <- 0 #makes all NAs into 0s
  Les[1,] <- c(Fert)
  class(Les) <- "leslie.matrix"
  Lambda <- popbio::eigen.analysis(Les)$lambda1
  L[i, 1:3] <- c(paste0(N2_wide$Treatment[i]), paste0(N2_wide$ID[i]), Lambda)
  
}


colnames(L)<-c("Treatment", "ID", "Lambda")

Data<-as.data.frame(L)

Data$Lambda<-as.numeric(as.character(Data$Lambda))

Lambda_N2 <-
  Data %>%
  dabest(Treatment, Lambda,
         idx = list(c("ev", "daf")))

## Lambda dabestr plot
Lambda_plot_N2 <- plot(mean_diff(Lambda_N2),
     #color.column = Treatment,
     palette = c('darkorange2','deepskyblue4'), axes.title.fontsize=13, 
                 rawplot.ylabel="Fitness", effsize.ylabel = '', rawplot.markersize=2, tick.fontsize = 9)

Lambda_plot_N2

### LRS

totalrep<-na.omit(as.data.frame.table(tapply(N2_all$value,list(N2_all$Treatment, N2_all$ID_ind),sum)))


names(totalrep)<-c("Treatment", "Replicate", "Totrep")

Totrep <-
  totalrep %>%
  dabest(Treatment, Totrep,
         idx = c('ev','daf'))

LRS <- plot(mean_diff(Totrep), palette = c('darkorange2','deepskyblue4'),
          #color.column = Treatment, #this remove the legend - make sure colours are correct though
          axes.title.fontsize=13, 
          rawplot.ylabel="Total reproduction", effsize.ylabel = "Mean difference", rawplot.markersize=2, tick.fontsize=9)
LRS
```