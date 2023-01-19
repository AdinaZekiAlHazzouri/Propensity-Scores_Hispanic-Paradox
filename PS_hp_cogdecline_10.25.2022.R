# Project: "PS Hispanic Paradox" manuscript 
# Last update: 10/25/2022

#install.packages(c("tidyverse","tidyr","nnet","nlme","lme4","optimx","ggeffects"))  
library(tidyverse)
library(haven)
library(tidyr)
library(nnet)
library(nlme)
library(lme4)
library(optimx)
library(ggeffects)
library(cowplot)

########################################################################################################################################################################################################
########################################################################################################################################################################################################

# Load dataset (merging PITCH and HRS)

data = read_csv("/Users/klk2131/Projects/PS_HP/Datasets/Merge_hrs_pitch.csv",
                col_types = cols(.default = col_double())) 
#rand <- read_dta("/Users/klk2131/Projects/PS_HP/Datasets/randhrs1992_2018v1.dta")
rand92<-rand[,c("hhidpn","r1wtresp")] #get baseline weights
names(rand92) <- toupper(names(rand92))
data<-left_join(data,rand92,by="HHIDPN") #add baseline weights to dataset

########################################################################################################################################################################################################
########################################################################################################################################################################################################

#(1) DATA PREPARATION: data cleaning and converting dataset from wide to long format

# 1. Clean up variables and run descriptives  

## Father and mother's education - create a missing data category (=5)
data = data %>% 
  mutate(rafeduc_r = ifelse(is.na(rafeduc_m) == TRUE, 5, rafeduc_m),
         rameduc_r = ifelse(is.na(rameduc_m) == TRUE, 5, rameduc_m),
         rafeduc_r = factor(rafeduc_r),
         rameduc_r = factor(rameduc_r))

## Age first smoking
  # Scale and center
  # Recode NA to 0
  # Indicator variable

  data = data %>% 
    mutate(agefstsmk_r = scale(agefstsmk), 
           agefstsmk_r = ifelse(is.na(agefstsmk_r), 0, agefstsmk_r), 
           never_smk = ifelse(is.na(agefstsmk), 1, 0)) 

      #for later to have original variables for table 1
      data$agefstsmk_desc<-data$agefstsmk

## Age first marriage
  # Scale and center
  # Recode NA to 0
  # Indicator variable

  data = data %>% 
    mutate(agefstm_r = scale(agefstm),
         agefstm_r = ifelse(is.na(agefstm_r) == TRUE, 0, agefstm_r),
         never_m = ifelse(is.na(agefstm), 1, 0))

    #for later to have original variables for table 1
    data$agefstm_desc<-data$agefstm

## Age first job
  # Scale and center
  # Recode NA to 0
  # Indicator variable
  
  data = data %>% 
    mutate(agefstjob_r = scale(agefstjob),
           agefstjob_r = ifelse(is.na(agefstjob_r) == TRUE & RAJJOBS == 0, 0, agefstjob_r),
           never_job = ifelse(is.na(agefstjob), 1, 0))
  
    #for later to have original variables for table 1
    data$agefstjob_desc<-data$agefstjob

## Birth year (birth cohort) 
#rescale birth year for later analyses
data = data %>% 
    mutate(RABYEAR_r = scale(RABYEAR, 1934, 10))

#Create complete case dataset for PS analyses
  data = data %>% filter(!is.na(RABYEAR_r),!is.na(GENDER), !is.na(rafeduc_r), !is.na(rameduc_r), !is.na(RAEDYRS), !is.na(height), !is.na(agefstjob_r), !is.na(agefstm_r), !is.na(agefstsmk_r))

  #n=22584
    
#for table 1
  prop.table(table(data$GENDER,data$group),margin=2)*100
  prop.table(table(data$rafeduc_r,data$group),margin=2)*100
  prop.table(table(data$rameduc_r,data$group),margin=2)*100
  prop.table(table(data$never_m,data$group),margin=2)*100
  prop.table(table(data$never_smk,data$group),margin=2)*100
  prop.table(table(data$never_job,data$group),margin=2)*100
  
  aggregate(data$RAEDYRS, list(data$group), mean)
  aggregate(data$RAEDYRS, list(data$group), sd)
  
  aggregate(data$RABYEAR, list(data$group), mean)
  aggregate(data$RABYEAR, list(data$group), sd)
  
  aggregate(data$height, list(data$group), mean)
  aggregate(data$height, list(data$group), sd)
  
  aggregate(data$agefstm_desc[!is.na(data$agefstm_desc)], list(data$group[!is.na(data$agefstm_desc)]), mean)
  aggregate(data$agefstm_desc[!is.na(data$agefstm_desc)], list(data$group[!is.na(data$agefstm_desc)]), sd)
  
  aggregate(data$agefstsmk_desc[!is.na(data$agefstsmk_desc)], list(data$group[!is.na(data$agefstsmk_desc)]), mean)
  aggregate(data$agefstsmk_desc[!is.na(data$agefstsmk_desc)], list(data$group[!is.na(data$agefstsmk_desc)]), sd)
  
  aggregate(data$agefstjob_desc[!is.na(data$agefstjob_desc)], list(data$group[!is.na(data$agefstjob_desc)]), mean)
  aggregate(data$agefstjob_desc[!is.na(data$agefstjob_desc)], list(data$group[!is.na(data$agefstjob_desc)]), sd)

########################################################################################################################################################################################################
########################################################################################################################################################################################################
    
#(2). FIT PROPENSITY SCORES  

## Set reference as US Born Non-Hispanic white; dp=PS denominator
data_ps = data %>% mutate(group = factor(group, levels = c("1", "2", "3"), 
                                         labels = c("dp_nhw", 
                                                    "dp_usbmxam", 
                                                    "dp_fbmxam")))

  table(data_ps$group)

## Multinomial regression to create propensity score model
fit1 = multinom(group ~ factor(GENDER) + RAEDYRS + factor(rafeduc_r) + factor(rameduc_r) +
                  height + agefstm_r + factor(never_m) + agefstsmk_r + factor(never_smk) +  
                  agefstjob_r + factor(never_job) + RABYEAR_r, data = data_ps)

## Calculate propensity scores
ps = predict(fit1, type = "probs")
  #check that sum of the os gives you 3x the sample size
    sum(table(ps)) #67752 = 3x the sample size

## Add dp to the dataset
ps = predict(fit1, type = "probs") %>% as.data.frame()
    #data checks
      ps$sum<-ps$dp_nhw+ps$dp_usbmxam+ps$dp_fbmxam
      summary(ps$sum) #data check: good, all probabilities for PS groups in each individual sum to 1
      sum(ps$dp_nhw)+sum(ps$dp_fbmxam)+sum(ps$dp_usbmxam) #data check: summing all the weights (value of the weights) for each group gives the original sample size, n=22584; value of each person's weights for each of the three exposure groups will sum to 1
      sum(table(ps$dp_nhw))+sum(table(ps$dp_usbmxam))+sum(table(ps$dp_fbmxam)) #summing all the people/# of observations given a weight for each group gives 3x the sample size (so each person given a weight for each group, 3x sample size because each person used 3 times to represent themselves under each exposure condition)

## Join propensity scores to original data
  data = cbind(data, ps)

  #distribution of PS for each group
  summary(round(data$dp_nhw, 2)) 
    hist(data$dp_nhw)
  summary(round(data$dp_fbmxam, 2))
    hist(data$dp_fbmxam)
  summary(round(data$dp_usbmxam, 2))
    hist(data$dp_usbmxam)
  
########################################################################################################################################################################################################
########################################################################################################################################################################################################
    
#(3). CREATE PS WEIGHTS

data = data %>% 
  mutate(w_hiscatg = case_when(group == 1 ~ 1/dp_nhw,
                               group == 2 ~ 1/dp_usbmxam,
                               group == 3 ~ 1/dp_fbmxam))

  #distribution of weights for total sample and each exposure group
  summary(round(data$w_hiscatg, 2)) #average weight (not stabilized) is around 5
  summary(round(data$w_hiscatg[data$group==1], 2)) #mean weight is 1.115, max is 333.80
  summary(round(data$w_hiscatg[data$group==2], 2)) #mean weight is 23.386, max is 1000.350
  summary(round(data$w_hiscatg[data$group==3], 2)) #mean weight is 70.68, max is 21663.11
  
#so these weights are unstable...need to  stabilize
  
########################################################################################################################################################################################################
########################################################################################################################################################################################################
  
#(4). CREATE STABILIZED WEIGHTS

#numerator weight approach 
  
# Set reference as US Born Non-Hispanic
data_np = data %>% mutate(group = factor(group, levels = c("1", "2", "3"), 
                                         labels = c("np_nhw", 
                                                    "np_usbmxam", 
                                                    "np_fbmxam")))

# Multinomial regression
fit2 = multinom(group ~ factor(GENDER) + RABYEAR_r + factor(rafeduc_r) + factor(rameduc_r), 
                data = data_np)

# Calculate propensity scores
np = predict(fit2, type = "probs") 

# Join propensity scores to original data
data = cbind(data, np)

# Calculate stabilized weights
data = data %>% 
  mutate(stblz.wts = case_when(group == 1 ~ np_nhw/dp_nhw,
                                group == 2 ~ np_usbmxam/dp_usbmxam,
                                group == 3 ~ np_fbmxam/dp_fbmxam))

summary(data$stblz.wts) #with weight stabilization mean of weights in the sample is now around 1 (1.02660) - which is good
quantile(data$stblz.wts, c(.01, .05, .25, .75, .95, .99))

#Supp Table 1: distribution of weights for total sample and each exposure group
  summary(round(data$stblz.wts, 2)) #stabilized weights have mean around 1.03 (when adding covariates to numerator)
  summary(data$stblz.wts)
  summary(round(data$stblz.wts[data$group==1], 2)) #mean is 1.01 max is 158.7
  summary(round(data$stblz.wts[data$group==2], 2)) #mean is 1.09 max is 61.1
  summary(round(data$stblz.wts[data$group==3], 2)) #mean is 1.25 max is 108.1

###########################################################################################################################################################
###########################################################################################################################################################

#(5). PLOTTING THE PROPENSITY SCORE(S)

#plot distribution of nh_white PS 

data %>% 
  #filter(dp_nhw<0.3) %>%
  #filter(dp_nhw>=0.30&dp_nhw<=0.70) %>%
  #filter(dp_nhw>0.7) %>%
  mutate(group = factor(group, levels = c("1", "2", "3"), 
                        labels = c("USBORN, Non-Hispanic white", 
                                   "USBORN, Mexican American", 
                                   "Foreign Born, Mexican American"))) %>% 
  ggplot(aes(x = dp_nhw)) + 
  geom_density(aes(fill = group), alpha = 0.8, show.legend=FALSE) + 
  labs(#title = "Propensity Score for being USBORN, Non-Hispanic White",
    x = "Propensity Score",
    y = "Density",
    fill = "Racial/Ethnic Group") + 
  viridis::scale_fill_viridis(discrete = T, direction = -1) + 
  theme_bw()

  # - - - - 
    #for suppl tables/figures using dp_fbmxam and dp_usbmxam use code above with filters below and replace x=dp_nhw with other ps as appropriate
    
  #plots for dist PS for being in usbmxam group
  data %>% 
    #filter(dp_usbmxam<0.3) %>%
    #filter(dp_usbmxam>=0.30&dp_usbmxam<0.7) %>%
    #filter(dp_usbmxam>0.7) %>%
    mutate(group = factor(group, levels = c("1", "2", "3"), 
                            labels = c("USBORN, Non-Hispanic white", 
                                       "USBORN, Mexican American", 
                                       "Foreign Born, Mexican American"))) %>% 
      ggplot(aes(x = dp_usbmxam)) + 
      geom_density(aes(fill = group), alpha = 0.8, show.legend=FALSE) + 
      labs(#title = "Propensity Score for being USBORN, Non-Hispanic White",
        x = "Propensity Score",
        y = "Density",
        fill = "Racial/Ethnic Group") + 
      viridis::scale_fill_viridis(discrete = T, direction = -1) + 
      theme_bw()  
    
    #plots for dist PS for being in fbmxam group
      data %>% 
        #filter(dp_fbmxam<0.3) %>%
        #filter(dp_fbmxam>=0.30&dp_fbmxam<=0.70) %>%
        #filter(dp_fbmxam>0.7) %>%
        mutate(group = factor(group, levels = c("1", "2", "3"), 
                              labels = c("USBORN, Non-Hispanic white", 
                                         "USBORN, Mexican American", 
                                         "Foreign Born, Mexican American"))) %>% 
        ggplot(aes(x = dp_nhw)) + 
        geom_density(aes(fill = group), alpha = 0.8, show.legend=FALSE) + 
        labs(#title = "Propensity Score for being USBORN, Non-Hispanic White",
          x = "Propensity Score",
          y = "Density",
          fill = "Racial/Ethnic Group") + 
        viridis::scale_fill_viridis(discrete = T, direction = -1) + 
        theme_bw()  
    # - - - - 
    
# Create indicators for group 2 (usbmxam) and group 3 (fbmxam)  
data = data %>% 
  mutate(usbmxam = ifelse(group == 2, 1, 0),
         fbmxam = ifelse(group == 3, 1, 0))

#PS Descriptives in text (Discussion)
table(data$group[data$dp_nhw>0.95]) #Results paragraph 2
table(data$group[data$dp_nhw>0.7]) #Discussion paragraph 3 
sum(table(data$group[data$dp_nhw>0.7])) #Discussion paragraph 3 
    
###########################################################################################################################################################
###########################################################################################################################################################

#(6). CREATE LONG FORM OF DATA

data_long = data %>%
  select(HHIDPN, R1WTRESP, group, GENDER, RAEDYRS, RABYEAR, RABYEAR_r, HACOHORT, rafeduc_r, rameduc_r, height, 
         agefstm_r, never_m, agefstsmk_r, never_smk, agefstjob_r, never_job, 
         dp_nhw:dp_fbmxam, w_hiscatg, np_nhw:np_fbmxam, stblz.wts, agefstjob_desc, agefstsmk_desc, agefstm_desc,
         usbmxam, fbmxam, firstAGE, firstIWYEAR,
         #A_g, F_g, G_g, H_g, J_g, K_g, L_g, M_g, N_g, O_g, P_g, Q_g, 
         AAGE:OAGE, AIWYEAR:OIWYEAR) %>% 
  pivot_longer(-c(HHIDPN:firstIWYEAR),
               names_to = c("wave", ".value"), 
               names_pattern = "(^[A-Q])(AGE|IWYEAR)") %>% 
  filter(!wave %in% c("B", "D", "C", "E")) %>%
  filter(!is.na(wave)) %>%
  filter(IWYEAR!=9998)
#n=131806

data_long_cog = data %>%
  select(HHIDPN,group, R1WTRESP,GENDER, RAEDYRS, RABYEAR, RABYEAR_r, HACOHORT, rafeduc_r, rameduc_r, height, 
         agefstm_r, never_m, agefstsmk_r, never_smk, agefstjob_r, never_job, 
         dp_nhw:dp_fbmxam, w_hiscatg, np_nhw:np_fbmxam, stblz.wts, agefstjob_desc, agefstsmk_desc, agefstm_desc,
         usbmxam, fbmxam, firstAGE, firstIWYEAR,  
         A_g, F_g, G_g, H_g, J_g, K_g, L_g, M_g, N_g, O_g, P_g, Q_g) %>% 
  pivot_longer(-c(HHIDPN:firstIWYEAR),
                names_to = c("PITCHwave", ".value"), 
                names_pattern = "(^[A-Q])(_g)") %>% 
  filter(!PITCHwave %in% c("B", "D", "C", "E")) %>%
  filter(!is.na(PITCHwave)) %>%
  rename(cog = "_g") %>% 
  # If wave P or Q, interview years are 2016 and 2018 respectively
  # Fill in missing interview years
  mutate(IWYEAR = case_when(PITCHwave == "P" ~ 2016,
                            PITCHwave == "Q" ~ 2018,
                            PITCHwave == "O" ~ 2014,
                            PITCHwave == "N" ~ 2012,
                            PITCHwave == "M" ~ 2010,
                            PITCHwave == "L" ~ 2008,
                            PITCHwave == "K" ~ 2006,
                            PITCHwave == "J" ~ 2004,
                            PITCHwave == "H" ~ 2002,
                            PITCHwave == "G" ~ 2000,
                            PITCHwave == "F" ~ 1998,
                            PITCHwave == "A" ~ 1994)) 


##Conversion from wave numbers to letters
#waveLookup_hrs <- list(
#  "HRS_W2"="A",
#  "HRS_W4"="F",
#  "HRS_W5"="G",
#  "HRS_W6"="H",
#  "HRS_W7"="J",
#  "HRS_W8"="K",
#  "HRS_W9"="L",
#  "HRS_W10"="M",
#  "HRS_W11"="N",
#  "HRS_W12"="O",
#  "HRS_W13"="P",
#  "HRS_W14"="Q"
#)


#could do full join and then when restrict to having cognitive scores this equals N if did a right instead (could do a right join because have to have cognition value (NA cog value will get booted from analysis/won't contribute))
data_long_test<-full_join(data_long, data_long_cog, by=c("HHIDPN","IWYEAR","group", "GENDER", "RAEDYRS", "RABYEAR", "RABYEAR_r", "HACOHORT", "rafeduc_r", "rameduc_r", "height", 
         "agefstm_r", "never_m", "agefstsmk_r", "never_smk", "agefstjob_r", "never_job", 
         "dp_nhw","dp_usbmxam","dp_fbmxam", "w_hiscatg", "np_nhw","np_usbmxam","np_fbmxam", "stblz.wts", "agefstjob_desc", "agefstsmk_desc", "agefstm_desc",
         "usbmxam", "fbmxam", "firstAGE", "firstIWYEAR","R1WTRESP"))

round(aggregate(data_long_test$cog, by=list(c(data_long_test$IWYEAR)), FUN="summary"),5) #cognitive tests only in waves where should be (1994, 1998, 2000-2018 by 2yr intervals)
  
## Create practice effect indicator for first PITCH score
data_long = data_long_test %>%
  filter(!is.na(cog)) %>% 
  group_by(HHIDPN) %>% 
  mutate(n = row_number(),
         first_cog = ifelse(n == 1, 1, 0)) %>% 
  select(-n) %>% 
  ungroup()

#n=124683

table(data_long$first_cog, data_long$IWYEAR)
round(aggregate(data_long$cog, by=list(c(data_long$IWYEAR)), FUN="summary"),5) #only have cognitive tests with years that make sense (no longer get cognitive tests with year=1992)

###########################################################################################################################################################
###########################################################################################################################################################

#(7). CLEAN AND SCALE AGE VARIABLE

table(data_long$AGE)
summary(data_long$AGE)

data_long = data_long %>% 
  group_by(HHIDPN) %>% 
  mutate(min_AGE = min(AGE, na.rm = T),
         min_IWYEAR = min(IWYEAR, na.rm = T),
         has_cog = ifelse(!is.na(cog),1,0),
         num_cog = sum(has_cog),
         age_test = case_when(is.na(AGE) & !is.na(cog) ~ firstAGE + (IWYEAR - firstIWYEAR),
                           TRUE ~ AGE),
         age_test_r = ((age_test-70)/10),
         AGE2 = case_when(is.na(AGE) & !is.na(cog) ~ firstAGE + (IWYEAR - firstIWYEAR),
                           TRUE ~ AGE),
         AGE2 = as.numeric(AGE2)) %>% 
  ungroup() %>% 
  mutate(AGE_R = scale(AGE2, center=70, scale = 10)) 
  #age is in units of decades and centered 

  summary(data_long$AGE_R)
  summary(data_long$age_test) #confirming that this centers age at 70 and divides it by 10 to scale to decades
  table(data_long$num_cog) #min time points contributed = 1 and maximum = 11 
  summary(data_long$num_cog) #median number time points contributed cognitive assessment 8

#For Table 1
  aggregate(data_long$num_cog, by=list(c(data_long$group)), FUN="median", na.rm=TRUE)
  aggregate(data_long$num_cog, by=list(c(data_long$group)), FUN="min", na.rm=TRUE)
  aggregate(data_long$num_cog, by=list(c(data_long$group)), FUN="max", na.rm=TRUE)
  summary(data_long$num_cog[data_long$group==1])
  summary(data_long$num_cog[data_long$group==2])
  summary(data_long$num_cog[data_long$group==3])
  

################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################
  
#(8). SUPPLEMENTAL DESCRIPTIVE TABLES
  
#supplemental descriptive table (Table S1)
#summary of weights
summary(data$stblz.wts)#total sample
summary(data$stblz.wts[data$group==1])#NHW
summary(data$stblz.wts[data$group==2])#USBMXAM
summary(data$stblz.wts[data$group==3])#FBMXAM

sum(table((data$dp_nhw[data$dp_nhw>0.05&data$dp_nhw<0.95]))) #5208
sum(table((data$group[data$dp_nhw>0.05&data$dp_nhw<0.95]))) #5208

#Total Sample
  #weights
  round(c(summary(data$stblz.wts), sd(data$stblz.wts)),2)#total sample
  round(quantile(data$stblz.wts, probs = c(0.01, 0.02, 0.05,0.10, 0.25,0.5, 0.75, 0.90, 0.95, 0.98, 0.99)),2)

#NHW Sample
  #weights
  round(c(summary(data$stblz.wts[data$group==1]), sd(data$stblz.wts[data$group==1])),2)#total sample
  round(quantile(data$stblz.wts[data$group==1], probs = c(0.01, 0.02, 0.05,0.10, 0.25,0.5, 0.75, 0.90, 0.95, 0.98, 0.99)),2)
  
#USB Sample
  #weights
  round(c(summary(data$stblz.wts[data$group==2]), sd(data$stblz.wts[data$group==2])),2)#total sample
  round(quantile(data$stblz.wts[data$group==2], probs = c(0.01, 0.02, 0.05,0.10, 0.25,0.5, 0.75, 0.90, 0.95, 0.98, 0.99)),2)
  
#FB Sample
  #weights
  round(c(summary(data$stblz.wts[data$group==3]), sd(data$stblz.wts[data$group==3])),2)#total sample
  round(quantile(data$stblz.wts[data$group==3], probs = c(0.01, 0.02, 0.05,0.10, 0.25,0.5, 0.75, 0.90, 0.95, 0.98, 0.99)),2)

########################################################################################################################################################################################################
########################################################################################################################################################################################################

#(9). REGRESSION MODELS IN TABLE/FIGURE 2 - ORIGINAL ANALYSIS (Total Sample)

## Model A: Unadjusted Model

fit0 = lmer(cog ~ factor(usbmxam) + factor(fbmxam) + AGE_R + 
                         AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + (AGE_R|HHIDPN) ,
                data = data_long)

fit0_coef = summary(fit0)$coefficients

#Results for Table 2
unadjusted<-rbind(fit0_coef[c(2,5,3,6)], fit0_coef[c(2,5,3,6)]-(1.96*fit0_coef[c(2,5,3,6),2]), fit0_coef[c(2,5,3,6)]+(1.96*fit0_coef[c(2,5,3,6),2]))
  rownames(unadjusted)<-c("beta", "l95ci", "u95ci")
  round(unadjusted,2)


  #Plotting results for Figure 1 (total sample)
  # Create data for visualizing results
  newdat = expand.grid(age = seq(50, 100, by = 1),
                       group = c(1:3)) %>% 
    mutate(AGE_R = scale(age, scale = 10),
           usbmxam = ifelse(group == 2, 1, 0),
           fbmxam = ifelse(group == 3, 1, 0),
           group = factor(group, levels = c(1:3), labels = c("US Born, Non-Hispanic White",
                                                             "US Born, Mexican American",
                                                             "Foreign Born, Mexican American")))
  
  newdat$cog = predict(fit0, newdat, re.form = ~0) #pull from total sample results
  
  # Plot Visualization
  unadj.plot<-newdat %>% 
    ggplot(aes(x = age, y = cog, color = group)) + 
    geom_line(size = 1, show.legend=FALSE) + 
    labs(x = "Age", 
         y = "Cognitive Z-Score",
         title = "Model 1: Unadjusted model") + 
    scale_y_continuous(limits=c(-2.5,2.5)) +
    viridis::scale_color_viridis(discrete = TRUE, direction = -1) +
    theme_bw()+
    theme(plot.title = element_text(size=11))
  unadj.plot

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## Model B: Conventional Adjustment Regression model

#models and results for Table 2
fit3 = lmer(cog ~ factor(rafeduc_r) + factor(rameduc_r) + factor(usbmxam) + factor(fbmxam) +  
              AGE_R + RAEDYRS + height + factor(GENDER) + RABYEAR_r + agefstm_r + 
              factor(never_m) + agefstsmk_r + factor(never_smk) + agefstjob_r + 
              factor(never_job) + factor(first_cog) + AGE_R * factor(usbmxam) + 
              AGE_R * factor(fbmxam) + (AGE_R|HHIDPN), 
            data = data_long)

fit3_coef = summary(fit3)$coefficients
convreg<-rbind(fit3_coef[c(10,24,11,25)], fit3_coef[c(10,24,11,25)]-(1.96*fit3_coef[c(10,24,11,25),2]), fit3_coef[c(10,24,11,25)]+(1.96*fit3_coef[c(10,24,11,25),2]))
  rownames(convreg)<-c("beta", "l95ci", "u95ci")
  round(convreg,2)
        
  #Plotting results for Figure 1 (total sample)
  # Create data for visualizing results
    newdat3 = expand.grid(age = seq(50, 100, by = 1),
                         group = c(1:3)) %>% 
      mutate(AGE_R = scale(age, scale = 10),
             usbmxam = ifelse(group == 2, 1, 0),
             fbmxam = ifelse(group == 3, 1, 0),
             group = factor(group, levels = c(1:3), labels = c("US Born, Non-Hispanic White",
                                                               "US Born, Mexican American",
                                                               "Foreign Born, Mexican American")))
    
    # Select covariates for predicted memory scores
    cov = data_long %>% select(HHIDPN, rafeduc_r, rameduc_r, RAEDYRS, height,
                               GENDER, RABYEAR_r, agefstm_r, never_m, agefstsmk_r,
                               never_smk, agefstjob_r, never_job, first_cog) %>% 
      filter(GENDER == 1, RABYEAR_r == 0, rafeduc_r == 4, rameduc_r == 4, height == 1.778, RAEDYRS==12, agefstm_r > -0.39 & agefstm_r < -0.019, first_cog == 0) #, rafeduc_r==1, rameduc_r==1, never_m==0, never_smk==0, never_job==0)
     
    newdat3 = bind_cols(newdat3, cov[1,])
  
    # Predict memory score
    newdat3$cog = predict(fit3, newdat3, re.form = ~0)
  
    # Plot Visualization
    fit3.t<-newdat3 %>% 
      ggplot(aes(x = age, y = cog, color = group)) + 
      geom_line(size = 1, show.legend=FALSE) + 
      labs(x = "Age", 
           y = "Cognitive Z-Score",
           title="Model 2: Conventional adjustment") + 
      scale_y_continuous(limits=c(-2.5,2.5)) +
      viridis::scale_color_viridis(discrete = TRUE, direction = -1) +
      theme_bw()+
      theme(plot.title = element_text(size=11))
    fit3.t

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## Model C: IPW with Stabilized Weights Model

reg.fit4 = lmer(cog ~ AGE_R + factor(usbmxam) + factor(fbmxam) + 
                  AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + 
                  RABYEAR_r + factor(GENDER) + factor(first_cog) + factor(rafeduc_r) + factor(rameduc_r) + (AGE_R|HHIDPN), 
                weights = stblz.wts, data = data_long)

  fit4_coef = summary(reg.fit4)$coefficients
  ipwreg<-rbind(fit4_coef[c(3,16,4,17)], fit4_coef[c(3,16,4,17)]-(1.96*fit4_coef[c(3,16,4,17),2]), fit4_coef[c(3,16,4,17)]+(1.96*fit4_coef[c(3,16,4,17),2]))
    rownames(ipwreg)<-c("beta", "l95ci", "u95ci")
    round(ipwreg,2)

    
  #Plotting results for Figure 1 (total sample)
  # Create data for visualizing results
    newdat4 = expand.grid(age = seq(50, 100, by = 1),
                        group = c(1:3)) %>% 
    mutate(AGE_R = scale(age, scale = 10),
           usbmxam = ifelse(group == 2, 1, 0),
           fbmxam = ifelse(group == 3, 1, 0),
           group = factor(group, levels = c(1:3), labels = c("US Born, Non-Hispanic white",
                                                             "US Born, Mexican American",
                                                             "Foreign Born, Mexican American")))
  # Select covariates for predicted memory scores
  cov = data_long %>% select(HHIDPN, GENDER, RABYEAR_r, rafeduc_r, rameduc_r, first_cog) %>% 
    filter(GENDER == 1, RABYEAR_r == 0, rafeduc_r == 4, rameduc_r == 4, first_cog == 0)
  
  #data_long$RABYEAR[data_long$RABYEAR_r==0] #predictions at birth year=1934
  newdat4 = bind_cols(newdat4, cov[1,])
  
  # Predict memory score
  newdat4$cog = predict(reg.fit4, newdat4, re.form = ~0)
  
  # Plot Visualization
  #First make a plot that extracts the legend (use this in the combined panel plot later)
  fit4.legend<-newdat4 %>% 
    ggplot(aes(x = age, y = cog, color = group)) + 
    geom_line(size = 1) + 
    labs(x = "Age", 
         y = "Cognitive Z-Score",
         group="Race/Ethnicity") + 
    viridis::scale_color_viridis(discrete = TRUE, direction = -1) +
    theme_bw()
  #make legend horizontal
    legend_b <- get_legend(
      fit4.legend + 
        guides(color = guide_legend(nrow = 1)) +
        theme(legend.title=element_blank(),
              legend.text = element_text(size=10),
              legend.position = c(1.5,.5)))
  #create plot for IPW model in the total (.t) sample
  fit4.t<-newdat4 %>% 
    ggplot(aes(x = age, y = cog, color = group)) + 
    geom_line(size = 1, show.legend=FALSE) + 
    labs(x = "Age", 
         y = "Cognitive Z-Score",
         title = "Model 3: IPTW model") + 
    scale_y_continuous(limits=c(-2.5,2.5)) +
    viridis::scale_color_viridis(discrete = TRUE, direction = -1) +
    theme_bw()+
    theme(plot.title = element_text(size=11))

#Save plots as a 3 panel plot for manuscript Figure 1
Fig1<-plot_grid(unadj.plot,fit3.t, fit4.t,legend_b,labels=c("A","B","C"),ncol=3,nrow=2, rel_heights=c(1,.1))
Fig1
#ggsave("/Users/klk2131/Projects/PS_HP/Figures/Figure1.png", width=9, height=4.5, dpi=500)


########################################################################################################################################################################################################
########################################################################################################################################################################################################

#(10). REGRESSION MODELS IN TABLE/FIGURE 2 - R&R ANALYSIS ; 3-WAY TRIMMING (CRUMP approach: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6395163/)

## Crump trimming
  data_long_crump<-data_long[data_long$dp_nhw>=0.05&data_long$dp_usbmxam>=0.05&data_long$dp_fbmxam>=0.05,]
      table(data_long_crump$group)
      table(data_long$group)
        summary(data_long_crump$dp_nhw)
        summary(data_long_crump$dp_usbmxam)
        summary(data_long_crump$dp_fbmxam)
      
      #looking at distribution of each PS post-trimming
        data_long_crump %>% 
        mutate(group = factor(group, levels = c("1", "2", "3"), 
                              labels = c("USBORN, Non-Hispanic white", 
                                             "USBORN, Mexican American", 
                                             "Foreign Born, Mexican American"))) %>% 
        ggplot(aes(x = dp_nhw)) + #replace dp_nhw with dp_usbmxam and dp_fbmxam to see dist of other ps
        geom_density(aes(fill = group), alpha = 0.8, show.legend=FALSE) + 
        labs(#title = "Propensity Score for being USBORN, Non-Hispanic White",
          x = "Propensity Score",
          y = "Density",
          fill = "Racial/Ethnic Group") + 
        viridis::scale_fill_viridis(discrete = T, direction = -1) + 
        theme_bw()
    
      #sample size/distributions post trimming
      table(data$group[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05])

      summary(data_long_crump$dp_nhw)
      summary(data_long_crump$dp_usbmxam)
      summary(data_long_crump$dp_fbmxam)
      
      summary(data_long$dp_nhw)
      summary(data_long$dp_usbmxam)
      summary(data_long$dp_fbmxam)
      
## Running analyses in trimmed sample

      ## Model A: Unadjusted Model 
      fit0.c = lmer(cog ~ factor(usbmxam) + factor(fbmxam) + AGE_R + 
                      AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + (AGE_R|HHIDPN) ,
                    data = data_long_crump)
      
      fit0_coef.c = summary(fit0.c)$coefficients

        #Results for Table 2
        unadjusted.c<-rbind(fit0_coef.c[c(2,5,3,6)], fit0_coef.c[c(2,5,3,6)]-(1.96*fit0_coef.c[c(2,5,3,6),2]), fit0_coef.c[c(2,5,3,6)]+(1.96*fit0_coef.c[c(2,5,3,6),2]))
        rownames(unadjusted.c)<-c("beta", "l95ci", "u95ci")
        round(unadjusted.c,2)

      ## Model B: Conventional Adjustment Regression model
      fit3.c = lmer(cog ~ factor(rafeduc_r) + factor(rameduc_r) + factor(usbmxam) + factor(fbmxam) +  
                      AGE_R + RAEDYRS + height + factor(GENDER) + RABYEAR_r + agefstm_r + 
                      factor(never_m) + agefstsmk_r + factor(never_smk) + agefstjob_r + 
                      factor(never_job) + factor(first_cog) + AGE_R * factor(usbmxam) + 
                      AGE_R * factor(fbmxam) + (AGE_R|HHIDPN), 
                    data = data_long_crump)
      
      fit3_coef.c = summary(fit3.c)$coefficients

        #Results for Table 2
        convreg.c<-rbind(fit3_coef.c[c(10,24,11,25)], fit3_coef.c[c(10,24,11,25)]-(1.96*fit3_coef.c[c(10,24,11,25),2]), fit3_coef.c[c(10,24,11,25)]+(1.96*fit3_coef.c[c(10,24,11,25),2]))
        rownames(convreg.c)<-c("beta", "l95ci", "u95ci")
        round(convreg.c,2)
      
      ## Model C: IPW with Stabilized Weights Model
      reg.fit4.c = lmer(cog ~ AGE_R + factor(usbmxam) + factor(fbmxam) + 
                          AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + 
                          RABYEAR_r + factor(GENDER) + factor(first_cog) + factor(rafeduc_r) + factor(rameduc_r) + (AGE_R|HHIDPN), 
                        weights = stblz.wts, data = data_long_crump)
      
      fit4_coef.c = summary(reg.fit4.c)$coefficients

        #Results for Table 2
        ipwreg.c<-rbind(fit4_coef.c[c(3,16,4,17)], fit4_coef.c[c(3,16,4,17)]-(1.96*fit4_coef.c[c(3,16,4,17),2]), fit4_coef.c[c(3,16,4,17)]+(1.96*fit4_coef.c[c(3,16,4,17),2]))
        rownames(ipwreg.c)<-c("beta", "l95ci", "u95ci")
        round(ipwreg.c,2)

########################################################################################################################################################################################################
########################################################################################################################################################################################################
      
#(11). SUPPLEMENTAL DESCRIPTIVES IN TRIMMED SAMPLE (SUPPLEMENTAL TABLE 2 & 4) 
      
## Supp. Table 4 (descriptives in Crump trimmed sample; total trimmed sample and by exposure group)
      #getting Ns of different strata of exposure/propensity score in the trimmed sample
       #RERUN/REPEAT for : GENDER, rafeduc_r, rameduc_r, rameduc_r, never_m, never_job, never_job
       #entire trimmed sample
       table(data$never_job[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==1]) #1=MALE, 2=FEMALE
       table(data$never_job[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==2]) #1=MALE, 2=FEMALE
       table(data$never_job[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==3]) #1=MALE, 2=FEMALE
       
       round(prop.table(table(data$never_job[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==1]))*100,0) #1=MALE, 2=FEMALE
       round(prop.table(table(data$never_job[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==2]))*100,0) #1=MALE, 2=FEMALE
       round(prop.table(table(data$never_job[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==3]))*100,0) #1=MALE, 2=FEMALE
       
      #Supp Table 4 education
      round(c(mean(data$RAEDYRS[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==1]),sd(data$RAEDYRS[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==1])),1)
      round(c(mean(data$RAEDYRS[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==2]),sd(data$RAEDYRS[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==2])),1)
      round(c(mean(data$RAEDYRS[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==3]),sd(data$RAEDYRS[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==3])),1)
      
      #birth year
      #Supp Table 4
      round(c(mean(data$RABYEAR[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==1]),sd(data$RABYEAR[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==1])),1)
      round(c(mean(data$RABYEAR[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==2]),sd(data$RABYEAR[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==2])),1)
      round(c(mean(data$RABYEAR[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==3]),sd(data$RABYEAR[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==3])),1)
      
      #height
      #Supp Table 4
      round(c(mean(data$height[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==1]),sd(data$height[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==1])),1)
      round(c(mean(data$height[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==2]),sd(data$height[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==2])),1)
      round(c(mean(data$height[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==3]),sd(data$height[data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==3])),1)
      
      #age first married
      #Supp Table 4
      round(c(mean(data$agefstm_desc[!is.na(data$agefstm_desc)&data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==1]),sd(data$agefstm_desc[!is.na(data$agefstm_desc)&data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==1])),1)
      round(c(mean(data$agefstm_desc[!is.na(data$agefstm_desc)&data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==2]),sd(data$agefstm_desc[!is.na(data$agefstm_desc)&data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==2])),1)
      round(c(mean(data$agefstm_desc[!is.na(data$agefstm_desc)&data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==3]),sd(data$agefstm_desc[!is.na(data$agefstm_desc)&data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==3])),1)

      #age first smoke
      #Supp Table 4
      round(c(mean(data$agefstsmk_desc[!is.na(data$agefstsmk_desc)&data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==1]),sd(data$agefstsmk_desc[!is.na(data$agefstsmk_desc)&data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==1])),1)
      round(c(mean(data$agefstsmk_desc[!is.na(data$agefstsmk_desc)&data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==2]),sd(data$agefstsmk_desc[!is.na(data$agefstsmk_desc)&data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==2])),1)
      round(c(mean(data$agefstsmk_desc[!is.na(data$agefstsmk_desc)&data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==3]),sd(data$agefstsmk_desc[!is.na(data$agefstsmk_desc)&data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==3])),1)

      #age first job
      #Supp Table 4
      round(c(mean(data$agefstjob_desc[!is.na(data$agefstjob_desc)&data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==1]),sd(data$agefstjob_desc[!is.na(data$agefstjob_desc)&data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==1])),1)
      round(c(mean(data$agefstjob_desc[!is.na(data$agefstjob_desc)&data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==2]),sd(data$agefstjob_desc[!is.na(data$agefstjob_desc)&data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==2])),1)
      round(c(mean(data$agefstjob_desc[!is.na(data$agefstjob_desc)&data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==3]),sd(data$agefstjob_desc[!is.na(data$agefstjob_desc)&data$dp_nhw>=0.05&data$dp_usbmxam>=0.05&data$dp_fbmxam>=0.05&data$group==3])),1)

      #median num cog
      summary(data_long$num_cog[data_long$group==1])
      summary(data_long$num_cog[data_long$group==2])
      summary(data_long$num_cog[data_long$group==3])
      summary(data_long_crump$num_cog)
      summary(data_long_crump$num_cog[data_long_crump$group==1])
      summary(data_long_crump$num_cog[data_long_crump$group==2])
      summary(data_long_crump$num_cog[data_long_crump$group==3])

########################################################################################################################################################################################################
########################################################################################################################################################################################################
      
#(12). SUPPLEMENTAL DESCRIPTIVES IN PROPENSITY SCORE STRATA (LOW,MED,HIGH): SUPPLEMENTAL TABLE 2
      #SAMPLE SIZES WITHIN STRATA OF OTHER PS (PS_FBMXAM, PS_USBMXAM): SUPPLEMENTAL TABLE 3
      
## Supplemental descriptive tables looking at PS by covariates and looking at balance in distrbution of covariates by exposure groups within propensity score strata  

  #(a) getting Ns of different strata of exposure/propensity score for Tables S2
      sum(table(unique(data_long$HHIDPN[data_long$group==1])))
      sum(table(unique(data_long$HHIDPN[data_long$group==2])))
      sum(table(unique(data_long$HHIDPN[data_long$group==3])))
      sum(table(unique(data_long$HHIDPN[data_long$group==1&data_long$dp_nhw<0.3])))
      sum(table(unique(data_long$HHIDPN[data_long$group==2&data_long$dp_nhw<0.3])))
      sum(table(unique(data_long$HHIDPN[data_long$group==3&data_long$dp_nhw<0.3])))
      sum(table(unique(data_long$HHIDPN[data_long$group==1&data_long$dp_nhw>=0.3&data_long$dp_nhw<=0.7])))
      sum(table(unique(data_long$HHIDPN[data_long$group==2&data_long$dp_nhw>=0.3&data_long$dp_nhw<=0.7])))
      sum(table(unique(data_long$HHIDPN[data_long$group==3&data_long$dp_nhw>=0.3&data_long$dp_nhw<=0.7])))
      sum(table(unique(data_long$HHIDPN[data_long$group==1&data_long$dp_nhw>0.7])))
      sum(table(unique(data_long$HHIDPN[data_long$group==2&data_long$dp_nhw>0.7])))
      sum(table(unique(data_long$HHIDPN[data_long$group==3&data_long$dp_nhw>0.7])))
      
      #Supp. Table 3 getting Ns of different strata of exposure/propensity score for supplemental tables (n within strata of PS fbmxam)
      sum(table(unique(data_long$HHIDPN[data_long$group==1])))
      sum(table(unique(data_long$HHIDPN[data_long$group==2])))
      sum(table(unique(data_long$HHIDPN[data_long$group==3])))
      sum(table(unique(data_long$HHIDPN[data_long$group==1&data_long$dp_fbmxam<0.3])))
      sum(table(unique(data_long$HHIDPN[data_long$group==2&data_long$dp_fbmxam<0.3])))
      sum(table(unique(data_long$HHIDPN[data_long$group==3&data_long$dp_fbmxam<0.3])))
      sum(table(unique(data_long$HHIDPN[data_long$group==1&data_long$dp_fbmxam>=0.3&data_long$dp_fbmxam<=0.7])))
      sum(table(unique(data_long$HHIDPN[data_long$group==2&data_long$dp_fbmxam>=0.3&data_long$dp_fbmxam<=0.7])))
      sum(table(unique(data_long$HHIDPN[data_long$group==3&data_long$dp_fbmxam>=0.3&data_long$dp_fbmxam<=0.7])))
      sum(table(unique(data_long$HHIDPN[data_long$group==1&data_long$dp_fbmxam>0.7])))
      sum(table(unique(data_long$HHIDPN[data_long$group==2&data_long$dp_fbmxam>0.7])))
      sum(table(unique(data_long$HHIDPN[data_long$group==3&data_long$dp_fbmxam>0.7])))
      
      #Supp Table 3 getting Ns of different strata of exposure/propensity score for supplemental tables (n within strata of PS usbmxam)
      sum(table(unique(data_long$HHIDPN[data_long$group==1])))
      sum(table(unique(data_long$HHIDPN[data_long$group==2])))
      sum(table(unique(data_long$HHIDPN[data_long$group==3])))
      sum(table(unique(data_long$HHIDPN[data_long$group==1&data_long$dp_usbmxam<0.3])))
      sum(table(unique(data_long$HHIDPN[data_long$group==2&data_long$dp_usbmxam<0.3])))
      sum(table(unique(data_long$HHIDPN[data_long$group==3&data_long$dp_usbmxam<0.3])))
      sum(table(unique(data_long$HHIDPN[data_long$group==1&data_long$dp_usbmxam>=0.3&data_long$dp_usbmxam<=0.7])))
      sum(table(unique(data_long$HHIDPN[data_long$group==2&data_long$dp_usbmxam>=0.3&data_long$dp_usbmxam<=0.7])))
      sum(table(unique(data_long$HHIDPN[data_long$group==3&data_long$dp_usbmxam>=0.3&data_long$dp_usbmxam<=0.7])))
      sum(table(unique(data_long$HHIDPN[data_long$group==1&data_long$dp_usbmxam>0.7])))
      sum(table(unique(data_long$HHIDPN[data_long$group==2&data_long$dp_usbmxam>0.7])))
      sum(table(unique(data_long$HHIDPN[data_long$group==3&data_long$dp_usbmxam>0.7])))
      
  #(b) characteristics within strata of PS in total sample (Supp Table 2)
      #Supp Table 2
      #getting Ns of different strata of exposure/propensity score for Tables S2 and S3
      #categorical predictors first
        #REPEAT ANALYSES FOR: GENDER, rafeduc_r, rameduc_r, rameduc_r, never_m, never_job, never_job

      #strata of PS for being NHW >0.05 to <0.3
      table(data$never_job[data$dp_nhw<0.3&data$group==1]) #1=MALE, 2=FEMALE
      table(data$never_job[data$dp_nhw<0.3&data$group==2]) #1=MALE, 2=FEMALE
      table(data$never_job[data$dp_nhw<0.3&data$group==3]) #1=MALE, 2=FEMALE
      
      round(prop.table(table(data$never_job[data$dp_nhw<0.3&data$group==1]))*100,0) #1=MALE, 2=FEMALE
      round(prop.table(table(data$never_job[data$dp_nhw<0.3&data$group==2]))*100,0) #1=MALE, 2=FEMALE
      round(prop.table(table(data$never_job[data$dp_nhw<0.3&data$group==3]))*100,0) #1=MALE, 2=FEMALE

      #strata of PS for being NHW >=0.3 to <=0.7
      table(data$never_job[data$dp_nhw>=0.3&data$dp_nhw<=0.7&data$group==1]) #1=MALE, 2=FEMALE
      table(data$never_job[data$dp_nhw>=0.3&data$dp_nhw<=0.7&data$group==2]) #1=MALE, 2=FEMALE
      table(data$never_job[data$dp_nhw>=0.3&data$dp_nhw<=0.7&data$group==3]) #1=MALE, 2=FEMALE
      
      round(prop.table(table(data$never_job[data$dp_nhw>=0.3&data$dp_nhw<=0.7&data$group==1]))*100,0) #1=MALE, 2=FEMALE
      round(prop.table(table(data$never_job[data$dp_nhw>=0.3&data$dp_nhw<=0.7&data$group==2]))*100,0) #1=MALE, 2=FEMALE
      round(prop.table(table(data$never_job[data$dp_nhw>=0.3&data$dp_nhw<=0.7&data$group==3]))*100,0) #1=MALE, 2=FEMALE
      
      #strata of PS for being NHW >0.7
      table(data$never_job[data$dp_nhw>0.7&data$group==1]) #1=MALE, 2=FEMALE
      table(data$never_job[data$dp_nhw>0.7&data$group==2]) #1=MALE, 2=FEMALE
      table(data$never_job[data$dp_nhw>0.7&data$group==3]) #1=MALE, 2=FEMALE
      
      round(prop.table(table(data$never_job[data$dp_nhw>0.7&data$group==1]))*100,0) #1=MALE, 2=FEMALE
      round(prop.table(table(data$never_job[data$dp_nhw>0.7&data$group==2]))*100,0) #1=MALE, 2=FEMALE
      round(prop.table(table(data$never_job[data$dp_nhw>0.7&data$group==3]))*100,0) #1=MALE, 2=FEMALE

    #Continuous covariates
      #For Supp Table 2 continuous covariates
      #create a PS categorical variable for descriptive analyses
      data$ps.strata<-ifelse(data$dp_nhw<0.3&data$group==1,"<0.3, NHW",
                             ifelse(data$dp_nhw>=0.3&data$dp_nhw<=0.7&data$group==1,"0.3-0.7, NHW",
                                    ifelse(data$dp_nhw>0.7&data$group==1,">0.7, NHW",
                                           ifelse(data$dp_nhw<0.3&data$group==2,"<0.3, USBMXAM",
                                                  ifelse(data$dp_nhw>=0.3&data$dp_nhw<=0.7&data$group==2,"0.3-0.7, USBMXAM",
                                                         ifelse(data$dp_nhw>0.7&data$group==2,">0.7, USBMXAM",
                                                                ifelse(data$dp_nhw<0.3&data$group==3,"<0.3, FBMXAM",
                                                                       ifelse(data$dp_nhw>=0.3&data$dp_nhw<=0.7&data$group==3,"0.3-0.7, FBMXAM",
                                                                              ifelse(data$dp_nhw>0.7&data$group==3,">0.7, FBMXAM",NA)))))))))
      data$ps.strata<-factor(data$ps.strata, levels = c("<0.3, NHW","<0.3, USBMXAM","<0.3, FBMXAM","0.3-0.7, NHW","0.3-0.7, USBMXAM","0.3-0.7, FBMXAM",">0.7, NHW",">0.7, USBMXAM",">0.7, FBMXAM"))
      table(data$ps.strata)
      
      #Yrs of education
      cbind(round(aggregate(data$RAEDYRS, by=list(c(data$ps.strata)), FUN="mean", na.rm=TRUE)[2],1), round(aggregate(data$RAEDYRS, by=list(c(data$ps.strata)), FUN="sd", na.rm=TRUE)[2],1))
      #birth year
      cbind(round(aggregate(data$RABYEAR, by=list(c(data$ps.strata)), FUN="mean", na.rm=TRUE)[2]), round(aggregate(data$RABYEAR, by=list(c(data$ps.strata)), FUN="sd", na.rm=TRUE)[2],1))
      #height
      cbind(round(aggregate(data$height, by=list(c(data$ps.strata)), FUN="mean", na.rm=TRUE)[2],1), round(aggregate(data$height, by=list(c(data$ps.strata)), FUN="sd", na.rm=TRUE)[2],1))
      #age first married
      cbind(round(aggregate(data$agefstm_desc, by=list(c(data$ps.strata)), FUN="mean", na.rm=TRUE)[2],1), round(aggregate(data$agefstm_desc, by=list(c(data$ps.strata)), FUN="sd", na.rm=TRUE)[2],1))
      #age first smoke
      cbind(round(aggregate(data$agefstsmk_desc, by=list(c(data$ps.strata)), FUN="mean", na.rm=TRUE)[2],1),round(aggregate(data$agefstsmk_desc, by=list(c(data$ps.strata)), FUN="sd", na.rm=TRUE)[2],1))
      #age first job
      cbind(round(aggregate(data$agefstjob_desc, by=list(c(data$ps.strata)), FUN="mean", na.rm=TRUE)[2],1),      round(aggregate(data$agefstjob_desc, by=list(c(data$ps.strata)), FUN="sd", na.rm=TRUE)[2],1))
      
      #median num cog
      summary(data_long$num_cog[data_long$group==1])
      summary(data_long$num_cog[data_long$group==2])
      summary(data_long$num_cog[data_long$group==3])
      summary(data_long$num_cog[data_long$dp_nhw<0.3&data_long$group==1])
      summary(data_long$num_cog[data_long$dp_nhw<0.3&data_long$group==2])
      summary(data_long$num_cog[data_long$dp_nhw<0.3&data_long$group==3])
      summary(data_long$num_cog[data_long$dp_nhw>=0.3&data_long$dp_nhw>=0.7&data_long$group==1])
      summary(data_long$num_cog[data_long$dp_nhw>=0.3&data_long$dp_nhw>=0.7&data_long$group==2])
      summary(data_long$num_cog[data_long$dp_nhw>=0.3&data_long$dp_nhw>=0.7&data_long$group==3])
      summary(data_long$num_cog[data_long$dp_nhw>0.7&data_long$group==1])
      summary(data_long$num_cog[data_long$dp_nhw>0.7&data_long$group==2])
      summary(data_long$num_cog[data_long$dp_nhw>0.7&data_long$group==3])
      
########################################################################################################################################################################################################
########################################################################################################################################################################################################

#(13) REGRESSION ANALYSES WITHIN STRATA OF THE PROPENSITY SCORE (TABLE 2)
      
## Unadjusted models
      #low PS for NHW
      fit0.lps= lmer(cog ~ factor(usbmxam) + factor(fbmxam) + AGE_R + 
                       AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + (AGE_R|HHIDPN) ,
                     data = subset(data_long, dp_nhw<0.3))
      
      #med PS for NHW
      fit0.mps= lmer(cog ~ factor(usbmxam) + factor(fbmxam) + AGE_R + 
                       AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + (AGE_R|HHIDPN) ,
                     data = subset(data_long, dp_nhw>=0.3&dp_nhw<=0.7))
      
      #high PS for NHW
      fit0.hps= lmer(cog ~ factor(usbmxam) + factor(fbmxam) + AGE_R + 
                       AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + (AGE_R|HHIDPN) ,
                     data = subset(data_long, dp_nhw>0.7))
      
      
      fit0_coef.l = summary(fit0.lps)$coefficients
      fit0_coef.m = summary(fit0.mps)$coefficients
      fit0_coef.h = summary(fit0.hps)$coefficients
      
      #Results for Table 2
      unadjusted.l<-rbind(fit0_coef.l[c(2,5,3,6)], fit0_coef.l[c(2,5,3,6)]-(1.96*fit0_coef.l[c(2,5,3,6),2]), fit0_coef.l[c(2,5,3,6)]+(1.96*fit0_coef.l[c(2,5,3,6),2]))
      rownames(unadjusted.l)<-c("beta", "l95ci", "u95ci")
      round(unadjusted.l,2)
      
      unadjusted.m<-rbind(fit0_coef.m[c(2,5,3,6)], fit0_coef.m[c(2,5,3,6)]-(1.96*fit0_coef.m[c(2,5,3,6),2]), fit0_coef.m[c(2,5,3,6)]+(1.96*fit0_coef.m[c(2,5,3,6),2]))
      rownames(unadjusted.m)<-c("beta", "l95ci", "u95ci")
      round(unadjusted.m,2)
      
      unadjusted.h<-rbind(fit0_coef.h[c(2,5,3,6)], fit0_coef.h[c(2,5,3,6)]-(1.96*fit0_coef.h[c(2,5,3,6),2]), fit0_coef.h[c(2,5,3,6)]+(1.96*fit0_coef.h[c(2,5,3,6),2]))
      rownames(unadjusted.h)<-c("beta", "l95ci", "u95ci")
      round(unadjusted.h,2)
      
## Conventional regression adjustment models
      #low PS for NHW
      fit3.lps = lmer(cog ~ factor(rafeduc_r) + factor(rameduc_r) + factor(usbmxam) + factor(fbmxam) +  
                        AGE_R + RAEDYRS + height + factor(GENDER) + RABYEAR_r + agefstm_r + 
                        factor(never_m) + agefstsmk_r + factor(never_smk) + agefstjob_r + 
                        factor(never_job) + factor(first_cog) + AGE_R * factor(usbmxam) + 
                        AGE_R * factor(fbmxam) + (AGE_R|HHIDPN), 
                      data = subset(data_long, dp_nhw<0.3))
      
      fit3_coef.l = summary(fit3.lps)$coefficients
      convreg.l<-rbind(fit3_coef.l[c(10,24,11,25)], fit3_coef.l[c(10,24,11,25)]-(1.96*fit3_coef.l[c(10,24,11,25),2]), fit3_coef.l[c(10,24,11,25)]+(1.96*fit3_coef.l[c(10,24,11,25),2]))
      rownames(convreg.l)<-c("beta", "l95ci", "u95ci")
      round(convreg.l,2)
      
      #med PS for NHW
      fit3.mps = lmer(cog ~ factor(rafeduc_r) + factor(rameduc_r) + factor(usbmxam) + factor(fbmxam) +  
                        AGE_R + RAEDYRS + height + factor(GENDER) + RABYEAR_r + agefstm_r + 
                        factor(never_m) + agefstsmk_r + factor(never_smk) + agefstjob_r + 
                        factor(never_job) + factor(first_cog) + AGE_R * factor(usbmxam) + 
                        AGE_R * factor(fbmxam) + (AGE_R|HHIDPN), 
                      data = subset(data_long, dp_nhw>=0.3&dp_nhw<=0.7))
      
      fit3_coef.m = summary(fit3.mps)$coefficients
      convreg.m<-rbind(fit3_coef.m[c(10,24,11,25)], fit3_coef.m[c(10,24,11,25)]-(1.96*fit3_coef.m[c(10,24,11,25),2]), fit3_coef.m[c(10,24,11,25)]+(1.96*fit3_coef.m[c(10,24,11,25),2]))
      rownames(convreg.m)<-c("beta", "l95ci", "u95ci")
      round(convreg.m,2)
      
      #high PS for NHW
      fit3.hps = lmer(cog ~ factor(rafeduc_r) + factor(rameduc_r) + factor(usbmxam) + factor(fbmxam) +  
                        AGE_R + RAEDYRS + height + factor(GENDER) + RABYEAR_r + agefstm_r + 
                        factor(never_m) + agefstsmk_r + factor(never_smk) + agefstjob_r + 
                        factor(never_job) + factor(first_cog) + AGE_R * factor(usbmxam) + 
                        AGE_R * factor(fbmxam) + (AGE_R|HHIDPN), 
                      data = subset(data_long, dp_nhw>0.7))
      
      fit3_coef.h = summary(fit3.hps)$coefficients
      convreg.h<-rbind(fit3_coef.h[c(10,24,11,25)], fit3_coef.h[c(10,24,11,25)]-(1.96*fit3_coef.h[c(10,24,11,25),2]), fit3_coef.h[c(10,24,11,25)]+(1.96*fit3_coef.h[c(10,24,11,25),2]))
      rownames(convreg.h)<-c("beta", "l95ci", "u95ci")
      round(convreg.h,2)
      
## IPTW models
      #low PS for NHW
      reg.fit4.lps = lmer(cog ~ AGE_R + factor(usbmxam) + factor(fbmxam) + 
                            AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + 
                            RABYEAR_r + factor(GENDER) + factor(first_cog) + factor(rafeduc_r) + factor(rameduc_r) + (AGE_R|HHIDPN), 
                          weights = stblz.wts, data = subset(data_long, dp_nhw<0.3))
      
      fit4_coef.l = summary(reg.fit4.lps)$coefficients
      ipwreg.l<-rbind(fit4_coef.l[c(3,16,4,17)], fit4_coef.l[c(3,16,4,17)]-(1.96*fit4_coef.l[c(3,16,4,17),2]), fit4_coef.l[c(3,16,4,17)]+(1.96*fit4_coef.l[c(3,16,4,17),2]))
      rownames(ipwreg.l)<-c("beta", "l95ci", "u95ci")
      round(ipwreg.l,2)
      
      #med PS for NHW
      reg.fit4.mps = lmer(cog ~ AGE_R + factor(usbmxam) + factor(fbmxam) + 
                            AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + 
                            RABYEAR_r + factor(GENDER) + factor(first_cog) + factor(rafeduc_r) + factor(rameduc_r) + (AGE_R|HHIDPN), 
                          weights = stblz.wts, data = subset(data_long, dp_nhw>=0.30&dp_nhw<=0.70))
      
      fit4_coef.m = summary(reg.fit4.mps)$coefficients
      ipwreg.m<-rbind(fit4_coef.m[c(3,16,4,17)], fit4_coef.m[c(3,16,4,17)]-(1.96*fit4_coef.m[c(3,16,4,17),2]), fit4_coef.m[c(3,16,4,17)]+(1.96*fit4_coef.m[c(3,16,4,17),2]))
      rownames(ipwreg.m)<-c("beta", "l95ci", "u95ci")
      round(ipwreg.m,2)
      
      #high PS for NHW
      reg.fit4.hps = lmer(cog ~ AGE_R + factor(usbmxam) + factor(fbmxam) + 
                            AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + 
                            RABYEAR_r + factor(GENDER) + factor(first_cog) + factor(rafeduc_r) + factor(rameduc_r) + (AGE_R|HHIDPN), 
                          weights = stblz.wts, data = subset(data_long, dp_nhw>0.70))
      
      fit4_coef.h = summary(reg.fit4.hps)$coefficients
      ipwreg.h<-rbind(fit4_coef.h[c(3,16,4,17)], fit4_coef.h[c(3,16,4,17)]-(1.96*fit4_coef.h[c(3,16,4,17),2]), fit4_coef.h[c(3,16,4,17)]+(1.96*fit4_coef.h[c(3,16,4,17),2]))
      rownames(ipwreg.h)<-c("beta", "l95ci", "u95ci")
      round(ipwreg.h,2)
      
########################################################################################################################################################################################################
########################################################################################################################################################################################################

#(14). SENSITIVITY ANALYSIS USING OTHER PS FOR STRATIFICATION: FBMXAM (Suppl Table S3)

    #SAMPLE SIZES
    table(data$group[data$dp_fbmxam&data$dp_fbmxam<0.3]) 
    table(data$group[data$dp_fbmxam>0.3&data$dp_fbmxam<0.7]) 
    table(data$group[data$dp_fbmxam>0.7])
    
## Model A: Unadjusted Models within strata of PS_FBMXAM
    
    #low PS fbmxam
    fit0.lps= lmer(cog ~ factor(usbmxam) + factor(fbmxam) + AGE_R + 
                     AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + (AGE_R|HHIDPN) ,
                   data = subset(data_long, dp_fbmxam<0.3))
    
    #med PS fbmxam
    fit0.mps= lmer(cog ~ factor(usbmxam) + factor(fbmxam) + AGE_R + 
                     AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + (AGE_R|HHIDPN) ,
                   data = subset(data_long, dp_fbmxam>=0.3&dp_fbmxam<=0.7))
    
    #high PS fbmxam
    fit0.hps= lmer(cog ~ factor(usbmxam) + factor(fbmxam) + AGE_R + 
                     AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + (AGE_R|HHIDPN) ,
                   data = subset(data_long, dp_fbmxam>0.7))
    
    
    fit0_coef.l = summary(fit0.lps)$coefficients
    fit0_coef.m = summary(fit0.mps)$coefficients
    fit0_coef.h = summary(fit0.hps)$coefficients
    
    #Results for Suppl Table 3
    unadjusted.l<-rbind(fit0_coef.l[c(2,5,3,6)], fit0_coef.l[c(2,5,3,6)]-(1.96*fit0_coef.l[c(2,5,3,6),2]), fit0_coef.l[c(2,5,3,6)]+(1.96*fit0_coef.l[c(2,5,3,6),2]))
    rownames(unadjusted.l)<-c("beta", "l95ci", "u95ci")
    round(unadjusted.l,2)
    
    unadjusted.m<-rbind(fit0_coef.m[c(2,5,3,6)], fit0_coef.m[c(2,5,3,6)]-(1.96*fit0_coef.m[c(2,5,3,6),2]), fit0_coef.m[c(2,5,3,6)]+(1.96*fit0_coef.m[c(2,5,3,6),2]))
    rownames(unadjusted.m)<-c("beta", "l95ci", "u95ci")
    round(unadjusted.m,2)
    
    unadjusted.h<-rbind(fit0_coef.h[c(2,5,3,6)], fit0_coef.h[c(2,5,3,6)]-(1.96*fit0_coef.h[c(2,5,3,6),2]), fit0_coef.h[c(2,5,3,6)]+(1.96*fit0_coef.h[c(2,5,3,6),2]))
    rownames(unadjusted.h)<-c("beta", "l95ci", "u95ci")
    round(unadjusted.h,2)
    
## Model B: Conventional Adjustment Regression models within strata of PS_FBMXAM
    
    #low PS fbmxam
    fit3.lps = lmer(cog ~ factor(rafeduc_r) + factor(rameduc_r) + factor(usbmxam) + factor(fbmxam) +  
                      AGE_R + RAEDYRS + height + factor(GENDER) + RABYEAR_r + agefstm_r + 
                      factor(never_m) + agefstsmk_r + factor(never_smk) + agefstjob_r + 
                      factor(never_job) + factor(first_cog) + AGE_R * factor(usbmxam) + 
                      AGE_R * factor(fbmxam) + (AGE_R|HHIDPN), 
                    data = subset(data_long, dp_fbmxam<0.3))
    
    fit3_coef.l = summary(fit3.lps)$coefficients
    convreg.l<-rbind(fit3_coef.l[c(10,24,11,25)], fit3_coef.l[c(10,24,11,25)]-(1.96*fit3_coef.l[c(10,24,11,25),2]), fit3_coef.l[c(10,24,11,25)]+(1.96*fit3_coef.l[c(10,24,11,25),2]))
    rownames(convreg.l)<-c("beta", "l95ci", "u95ci")
    round(convreg.l,2)
    
    #med PS fbmxam
    fit3.mps = lmer(cog ~ factor(rafeduc_r) + factor(rameduc_r) + factor(usbmxam) + factor(fbmxam) +  
                      AGE_R + RAEDYRS + height + factor(GENDER) + RABYEAR_r + agefstm_r + 
                      factor(never_m) + agefstsmk_r + factor(never_smk) + agefstjob_r + 
                      factor(never_job) + factor(first_cog) + AGE_R * factor(usbmxam) + 
                      AGE_R * factor(fbmxam) + (AGE_R|HHIDPN), 
                    data = subset(data_long, dp_fbmxam>=0.3&dp_fbmxam<=0.7))
    
    fit3_coef.m = summary(fit3.mps)$coefficients
    convreg.m<-rbind(fit3_coef.m[c(10,24,11,25)], fit3_coef.m[c(10,24,11,25)]-(1.96*fit3_coef.m[c(10,24,11,25),2]), fit3_coef.m[c(10,24,11,25)]+(1.96*fit3_coef.m[c(10,24,11,25),2]))
    rownames(convreg.m)<-c("beta", "l95ci", "u95ci")
    round(convreg.m,2)
    
    #high PS fbmxam
    fit3.hps = lmer(cog ~ factor(rafeduc_r) + factor(rameduc_r) + factor(usbmxam) + factor(fbmxam) +  
                      AGE_R + RAEDYRS + height + factor(GENDER) + RABYEAR_r + agefstm_r + 
                      factor(never_m) + agefstsmk_r + factor(never_smk) + agefstjob_r + 
                      factor(never_job) + factor(first_cog) + AGE_R * factor(usbmxam) + 
                      AGE_R * factor(fbmxam) + (AGE_R|HHIDPN), 
                    data = subset(data_long, dp_fbmxam>0.7))
    
    fit3_coef.h = summary(fit3.hps)$coefficients
    convreg.h<-rbind(fit3_coef.h[c(10,24,11,25)], fit3_coef.h[c(10,24,11,25)]-(1.96*fit3_coef.h[c(10,24,11,25),2]), fit3_coef.h[c(10,24,11,25)]+(1.96*fit3_coef.h[c(10,24,11,25),2]))
    rownames(convreg.h)<-c("beta", "l95ci", "u95ci")
    round(convreg.h,2)
    
## Model C: IPW with Stabilized Weights Models within strata of PS_FBMXAM
    
    #low PS fbmxam
    reg.fit4.lps = lmer(cog ~ AGE_R + factor(usbmxam) + factor(fbmxam) + 
                          AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + 
                          RABYEAR_r + factor(GENDER) + factor(first_cog) + factor(rafeduc_r) + factor(rameduc_r) + (AGE_R|HHIDPN), 
                        weights = stblz.wts, data = subset(data_long, dp_fbmxam<0.3))
    
    fit4_coef.l = summary(reg.fit4.lps)$coefficients
    ipwreg.l<-rbind(fit4_coef.l[c(3,16,4,17)], fit4_coef.l[c(3,16,4,17)]-(1.96*fit4_coef.l[c(3,16,4,17),2]), fit4_coef.l[c(3,16,4,17)]+(1.96*fit4_coef.l[c(3,16,4,17),2]))
    rownames(ipwreg.l)<-c("beta", "l95ci", "u95ci")
    round(ipwreg.l,2)
    
    #med PS fbmxam
    reg.fit4.mps = lmer(cog ~ AGE_R + factor(usbmxam) + factor(fbmxam) + 
                          AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + 
                          RABYEAR_r + factor(GENDER) + factor(first_cog) + factor(rafeduc_r) + factor(rameduc_r) + (AGE_R|HHIDPN), 
                        weights = stblz.wts, data = subset(data_long, dp_fbmxam>=0.30&dp_fbmxam<=0.70))
    
    fit4_coef.m = summary(reg.fit4.mps)$coefficients
    ipwreg.m<-rbind(fit4_coef.m[c(3,16,4,17)], fit4_coef.m[c(3,16,4,17)]-(1.96*fit4_coef.m[c(3,16,4,17),2]), fit4_coef.m[c(3,16,4,17)]+(1.96*fit4_coef.m[c(3,16,4,17),2]))
    rownames(ipwreg.m)<-c("beta", "l95ci", "u95ci")
    round(ipwreg.m,2)
    
    #high PS fbmxam
    reg.fit4.hps = lmer(cog ~ AGE_R + factor(usbmxam) + factor(fbmxam) + 
                          AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + 
                          RABYEAR_r + factor(GENDER) + factor(first_cog) + factor(rafeduc_r) + factor(rameduc_r) + (AGE_R|HHIDPN), 
                        weights = stblz.wts, data = subset(data_long, dp_fbmxam>0.70))
    
    fit4_coef.h = summary(reg.fit4.hps)$coefficients
    ipwreg.h<-rbind(fit4_coef.h[c(3,16,4,17)], fit4_coef.h[c(3,16,4,17)]-(1.96*fit4_coef.h[c(3,16,4,17),2]), fit4_coef.h[c(3,16,4,17)]+(1.96*fit4_coef.h[c(3,16,4,17),2]))
    rownames(ipwreg.h)<-c("beta", "l95ci", "u95ci")
    round(ipwreg.h,2)
    
########################################################################################################################################################################################################
########################################################################################################################################################################################################

#(15). SENSITIVITY ANALYSIS USING OTHER PS FOR STRATIFICATION: USBMXAM (Suppl Table S3)

    #SAMPLE SIZES
    table(data$group[data$dp_usbmxam<0.3]) 
    table(data$group[data$dp_usbmxam>0.3]) 
    table(data$group[data$dp_usbmxam>0.7]) #cannot run analyses here/n=1
    
## Model A: Unadjusted Model 

    #low PS usbmxam
    fit0.lps= lmer(cog ~ factor(usbmxam) + factor(fbmxam) + AGE_R + 
                     AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + (AGE_R|HHIDPN) ,
                   data = subset(data_long, dp_usbmxam<0.3))
    
    #med PS usbmxam
    fit0.mps= lmer(cog ~ factor(usbmxam) + factor(fbmxam) + AGE_R + 
                     AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + (AGE_R|HHIDPN) ,
                   data = subset(data_long, dp_usbmxam>=0.3&dp_usbmxam<=0.7))

    
    fit0_coef.l = summary(fit0.lps)$coefficients
    fit0_coef.m = summary(fit0.mps)$coefficients

    #Results for Supp Table  3
    unadjusted.l<-rbind(fit0_coef.l[c(2,5,3,6)], fit0_coef.l[c(2,5,3,6)]-(1.96*fit0_coef.l[c(2,5,3,6),2]), fit0_coef.l[c(2,5,3,6)]+(1.96*fit0_coef.l[c(2,5,3,6),2]))
    rownames(unadjusted.l)<-c("beta", "l95ci", "u95ci")
    round(unadjusted.l,2)
    
    unadjusted.m<-rbind(fit0_coef.m[c(2,5,3,6)], fit0_coef.m[c(2,5,3,6)]-(1.96*fit0_coef.m[c(2,5,3,6),2]), fit0_coef.m[c(2,5,3,6)]+(1.96*fit0_coef.m[c(2,5,3,6),2]))
    rownames(unadjusted.m)<-c("beta", "l95ci", "u95ci")
    round(unadjusted.m,2)
    
## Model B: Conventional Adjustment Regression model
    
    #low PS usbmxam
    fit3.lps = lmer(cog ~ factor(rafeduc_r) + factor(rameduc_r) + factor(usbmxam) + factor(fbmxam) +  
                      AGE_R + RAEDYRS + height + factor(GENDER) + RABYEAR_r + agefstm_r + 
                      factor(never_m) + agefstsmk_r + factor(never_smk) + agefstjob_r + 
                      factor(never_job) + factor(first_cog) + AGE_R * factor(usbmxam) + 
                      AGE_R * factor(fbmxam) + (AGE_R|HHIDPN), 
                    data = subset(data_long, dp_usbmxam<0.3))
    
    fit3_coef.l = summary(fit3.lps)$coefficients
    convreg.l<-rbind(fit3_coef.l[c(10,24,11,25)], fit3_coef.l[c(10,24,11,25)]-(1.96*fit3_coef.l[c(10,24,11,25),2]), fit3_coef.l[c(10,24,11,25)]+(1.96*fit3_coef.l[c(10,24,11,25),2]))
    rownames(convreg.l)<-c("beta", "l95ci", "u95ci")
    round(convreg.l,2)
    
    #med PS usbmxam
    fit3.mps = lmer(cog ~ factor(rafeduc_r) + factor(rameduc_r) + factor(usbmxam) + factor(fbmxam) +  
                      AGE_R + RAEDYRS + height + factor(GENDER) + RABYEAR_r + agefstm_r + 
                      factor(never_m) + agefstsmk_r + factor(never_smk) + agefstjob_r + 
                      factor(never_job) + factor(first_cog) + AGE_R * factor(usbmxam) + 
                      AGE_R * factor(fbmxam) + (AGE_R|HHIDPN), 
                    data = subset(data_long, dp_usbmxam>=0.3&dp_usbmxam<=0.7))
    
    fit3_coef.m = summary(fit3.mps)$coefficients
    convreg.m<-rbind(fit3_coef.m[c(10,24,11,25)], fit3_coef.m[c(10,24,11,25)]-(1.96*fit3_coef.m[c(10,24,11,25),2]), fit3_coef.m[c(10,24,11,25)]+(1.96*fit3_coef.m[c(10,24,11,25),2]))
    rownames(convreg.m)<-c("beta", "l95ci", "u95ci")
    round(convreg.m,2)
    
## Model C: IPW with Stabilized Weights Model

    #low PS usbmxam
    reg.fit4.lps = lmer(cog ~ AGE_R + factor(usbmxam) + factor(fbmxam) + 
                          AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + 
                          RABYEAR_r + factor(GENDER) + factor(first_cog) + factor(rafeduc_r) + factor(rameduc_r) + (AGE_R|HHIDPN), 
                        weights = stblz.wts, data = subset(data_long, dp_usbmxam<0.3))
    
    fit4_coef.l = summary(reg.fit4.lps)$coefficients
    ipwreg.l<-rbind(fit4_coef.l[c(3,16,4,17)], fit4_coef.l[c(3,16,4,17)]-(1.96*fit4_coef.l[c(3,16,4,17),2]), fit4_coef.l[c(3,16,4,17)]+(1.96*fit4_coef.l[c(3,16,4,17),2]))
    rownames(ipwreg.l)<-c("beta", "l95ci", "u95ci")
    round(ipwreg.l,2)
    
    #med PS usbmxam
    reg.fit4.mps = lmer(cog ~ AGE_R + factor(usbmxam) + factor(fbmxam) + 
                          AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + 
                          RABYEAR_r + factor(GENDER) + factor(first_cog) + factor(rafeduc_r) + factor(rameduc_r) + (AGE_R|HHIDPN), 
                        weights = stblz.wts, data = subset(data_long, dp_usbmxam>=0.30&dp_usbmxam<=0.70))
    
    fit4_coef.m = summary(reg.fit4.mps)$coefficients
    ipwreg.m<-rbind(fit4_coef.m[c(3,16,4,17)], fit4_coef.m[c(3,16,4,17)]-(1.96*fit4_coef.m[c(3,16,4,17),2]), fit4_coef.m[c(3,16,4,17)]+(1.96*fit4_coef.m[c(3,16,4,17),2]))
    rownames(ipwreg.m)<-c("beta", "l95ci", "u95ci")
    round(ipwreg.m,2)
    
########################################################################################################################################################################################################
########################################################################################################################################################################################################
#(16). MULTINOMIAL MATCHING (MATCH WEIGHTS FOR 3-WAY MATCHING)

#ref for approach: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5378668/
#ref: https://sejdemyr.github.io/r-tutorials/statistics/tutorial8.html
#https://cran.r-project.org/web/packages/MatchIt/vignettes/MatchIt.html
#With 1:1 nearest neighbor matching without replacement, excluding the matching weights does not change the estimates. For all other forms of matching, they are required, so we recommend always including them for consistency.
    
library(TriMatch)
library(MatchIt)
library(tableone)
library("cobalt") #for love plot to check covariate balance
library("sandwich")
library("lmtest")
library("VGAM")

    
#first run main analysis file up through creation of "data_ps"dataset and rename this dataset data_ps_matching; then start here
 #rerun lines 0-126 (Section 1 only)
 data_ps_matching<-data_ps
 
## Function to add generalized PS to dataset 
      #(creates PS for the three exposure groups as done above)
      #(later we can plot and see the PS created here as the same as dp_nhw, dp_usbmxam, and dp_fbmxam above)
    AddGPS <- function(data, formula, family = multinomial(), psPrefix = "PS_") {
      ## Fit multinomial logistic regression
      resVglm <- multinom(formula = group ~ factor(GENDER) + RAEDYRS + factor(rafeduc_r) + factor(rameduc_r) +
                            height + agefstm_r + factor(never_m) + agefstsmk_r + factor(never_smk) +  
                            agefstjob_r + factor(never_job) + RABYEAR_r, data = data_ps_matching, family = family)
      ## Calculate PS
      psData <- as.data.frame(predict(resVglm, type = "probs"))
      names(psData) <- paste0(psPrefix, names(psData))
      cbind(data, psData)
    }
    data_ps_matching <- AddGPS(data = data_ps_matching, # dataset
                               ## Propensity score model for multinomial regression
                               formula = group ~ factor(GENDER) + RAEDYRS + factor(rafeduc_r) + factor(rameduc_r) +
                                 height + agefstm_r + factor(never_m) + agefstsmk_r + factor(never_smk) +  
                                 agefstjob_r + factor(never_job) + RABYEAR_r)
    
    data_ps_matching$treat<-data_ps_matching$group
    table(data_ps_matching$treat)

    #2. CREATE WEIGHTS
    
    ## Function to add matching weight as mw to dataset
    AddMwToData <- function(data, txVar, txLevels, psPrefix = "PS_") {
      ## Treatment indicator data frame (any number of groups allowed)
      dfAssign <- as.data.frame(lapply(txLevels, function(tx_k) {
        as.numeric(data[txVar] == tx_k)
      }))
      ## Name of PS variables
      psVars <- paste0(psPrefix, txLevels)
      ## Pick denominator (PS for assigned treatment)
      data$PS_assign <- rowSums(data[psVars] * dfAssign)
      ## Pick numerator
      data$PS_min <- do.call(pmin, data[psVars])
      ## Calculate the IPTW
      data$iptw <- 1 / data$PS_assign
      ## Calculate the matching weight
      data$mw <- exp(log(data$PS_min) - log(data$PS_assign))
      ## Return the whole data
      data
    }
    
    ## Add IPTW and MW
    data_ps_matching<- AddMwToData(data = data_ps_matching, # dataset
                                    txVar = "treat", # treatment variable name
                                    txLevels = c("dp_nhw", "dp_usbmxam", "dp_fbmxam")) # treatment levels
    
    #All analyses afterward should be proceeded as weighted analyses. 
    
## Data checks to do:
    #(a) Plot distrbution of PS_dp_nw and make sure plot looks like dp_nhw from main analyses -- YES
    #(b) Plot distrbution of PS_dp_usbmxam and make sure plot looks like dp_usbmxam from main analyses -- YES
    #(c) Plot distrbution of PS_dp_fbmxam and make sure plot looks like dp_fbmxam from main analyses -- YES
      #creating a plot using e.g., PS_dp_nhw gives same plot as in main analyses, so creation of PS (above) is the same here
      #PS doesn't change, but weighting schematic does (to create a match weight vs an IPW) 
    
    #plot match weight and see how it looks
      data_ps_matching %>% 
      #filter(dp_nhw>0.05&dp_nhw<0.95) %>%
      #filter(dp_nhw>0.05&dp_nhw<0.3) %>%
      #filter(dp_nhw>=0.30&dp_nhw<=0.70) %>%
      #filter(dp_nhw>0.7&dp_nhw<0.95) %>%
      mutate(group = factor(group, levels = c("1", "2", "3"), 
                            labels = c("USBORN, Non-Hispanic white", 
                                       "USBORN, Mexican American", 
                                       "Foreign Born, Mexican American"))) %>% 
      ggplot(aes(x = mw)) + 
      geom_density(aes(fill = treat), alpha = 0.8, show.legend=FALSE) + 
      labs(#title = "Propensity Score for being USBORN, Non-Hispanic White",
        x = "Propensity Score",
        y = "Density",
        fill = "Racial/Ethnic Group") + 
      viridis::scale_fill_viridis(discrete = T, direction = -1) + 
      theme_bw()
    
## Proceed as before: create long dataset, clean up cognitive vars, create age variable  
    data_long_matching = data_ps_matching %>%
      select(HHIDPN, R1WTRESP, group, GENDER, RAEDYRS, RABYEAR, RABYEAR_r, HACOHORT, rafeduc_r, rameduc_r, height, 
             agefstm_r, never_m, agefstsmk_r, never_smk, agefstjob_r, never_job, 
             PS_dp_nhw, PS_dp_usbmxam, PS_dp_fbmxam, PS_assign, PS_min, iptw, mw, treat,firstAGE,firstIWYEAR,
             #A_g, F_g, G_g, H_g, J_g, K_g, L_g, M_g, N_g, O_g, P_g, Q_g, 
             AAGE:OAGE, AIWYEAR:OIWYEAR) %>% 
      pivot_longer(-c(HHIDPN:firstIWYEAR),
                   names_to = c("wave", ".value"), 
                   names_pattern = "(^[A-Q])(AGE|IWYEAR)") %>% 
      filter(!wave %in% c("B", "D", "C", "E")) %>%
      filter(!is.na(wave)) %>%
      filter(IWYEAR!=9998)
    #131806
    
    data_long_cog_matching = data_ps_matching %>%
      select(HHIDPN,group, R1WTRESP,GENDER, RAEDYRS, RABYEAR, RABYEAR_r, HACOHORT, rafeduc_r, rameduc_r, height, 
             agefstm_r, never_m, agefstsmk_r, never_smk, agefstjob_r, never_job, 
             PS_dp_nhw, PS_dp_usbmxam, PS_dp_fbmxam, PS_assign, PS_min, iptw, mw, treat,firstAGE,firstIWYEAR,
             A_g, F_g, G_g, H_g, J_g, K_g, L_g, M_g, N_g, O_g, P_g, Q_g) %>% 
      pivot_longer(-c(HHIDPN:firstIWYEAR),
                   names_to = c("PITCHwave", ".value"), 
                   names_pattern = "(^[A-Q])(_g)") %>% 
      filter(!PITCHwave %in% c("B", "D", "C", "E")) %>%
      filter(!is.na(PITCHwave)) %>%
      rename(cog = "_g") %>% 
      # If wave P or Q, interview years are 2016 and 2018 respectively
      # Fill in missing interview years
      mutate(IWYEAR = case_when(PITCHwave == "P" ~ 2016,
                                PITCHwave == "Q" ~ 2018,
                                PITCHwave == "O" ~ 2014,
                                PITCHwave == "N" ~ 2012,
                                PITCHwave == "M" ~ 2010,
                                PITCHwave == "L" ~ 2008,
                                PITCHwave == "K" ~ 2006,
                                PITCHwave == "J" ~ 2004,
                                PITCHwave == "H" ~ 2002,
                                PITCHwave == "G" ~ 2000,
                                PITCHwave == "F" ~ 1998,
                                PITCHwave == "A" ~ 1994)) 
    
    #could do full join and then when restrict to having cognitive scores this equals N if did a right instead (could do a right join because have to have cognition value (NA cog value will get booted from analysis/won't contribute))
    data_long_test<-full_join(data_long_matching, data_long_cog_matching, by=c("HHIDPN","IWYEAR","group", "GENDER", "RAEDYRS", "RABYEAR", "RABYEAR_r", "HACOHORT", "rafeduc_r", "rameduc_r", "height", 
                                                                               "agefstm_r", "never_m", "agefstsmk_r", "never_smk", "agefstjob_r", "never_job", 
                                                                               "PS_dp_nhw", "PS_dp_usbmxam", "PS_dp_fbmxam", "PS_assign", "PS_min", "iptw", "mw", "treat",
                                                                               "firstAGE","firstIWYEAR","R1WTRESP"))
    
    round(aggregate(data_long_test$cog, by=list(c(data_long_test$IWYEAR)), FUN="summary"),5) #cognitive tests only in waves where should be (1994, 1998, 2000-2018)
    table(data_long$IWYEAR)
    
    ## Create practice effect indicator for first PITCH score
    data_long_matching = data_long_test %>%
      filter(!is.na(cog)) %>% 
      group_by(HHIDPN) %>% 
      mutate(n = row_number(),
             first_cog = ifelse(n == 1, 1, 0)) %>% 
      select(-n) %>% 
      ungroup()
    
    #n=124683
    
    table(data_long$first_cog, data_long$IWYEAR)
    round(aggregate(data_long$cog, by=list(c(data_long$IWYEAR)), FUN="summary"),5) #only have cognitive tests with years that make sense (no longer get cognitive tests with year=1992)
    
    data_long_matching = data_long_matching %>% 
      group_by(HHIDPN) %>% 
      mutate(min_AGE = min(AGE, na.rm = T),
             min_IWYEAR = min(IWYEAR, na.rm = T),
             has_cog = ifelse(!is.na(cog),1,0),
             num_cog = sum(has_cog),
             age_test = case_when(is.na(AGE) & !is.na(cog) ~ firstAGE + (IWYEAR - firstIWYEAR),
                                  TRUE ~ AGE),
             age_test_r = ((age_test-70)/10),
             AGE2 = case_when(is.na(AGE) & !is.na(cog) ~ firstAGE + (IWYEAR - firstIWYEAR),
                              TRUE ~ AGE),
             AGE2 = as.numeric(AGE2)) %>% 
      ungroup() %>% 
      mutate(AGE_R = scale(AGE2, center=70, scale = 10)) 

## Run analysis weighted by matching weights
  #First descriptive/plots to check balance in match weighted sample
  #Then regression analyses
  
 #(A) plotting balance post-matching
    library(tableone)
    data_ps_matching$RABYEAR_r<-as.numeric(data_ps_matching$RABYEAR_r)
    data_ps_matching$agefstsmk_r<-as.numeric(data_ps_matching$agefstsmk_r)
    data_ps_matching$agefstm_r<-as.numeric(data_ps_matching$agefstm_r)
    data_ps_matching$agefstjob_r<-as.numeric(data_ps_matching$agefstjob_r)
    
    data_long_matching$RABYEAR_r<-as.numeric(data_long_matching$RABYEAR_r)
    data_long_matching$agefstsmk_r<-as.numeric(data_long_matching$agefstsmk_r)
    data_long_matching$agefstm_r<-as.numeric(data_long_matching$agefstm_r)
    data_long_matching$agefstjob_r<-as.numeric(data_long_matching$agefstjob_r)
    
    covariates <- c("GENDER", "RAEDYRS", "RABYEAR_r", "rafeduc_r",
                    "rameduc_r", "height", "agefstm_r", "agefstsmk_r",
                    "agefstjob_r", "never_m", "never_smk","never_job")
    tab1Unadj <- CreateTableOne(vars = covariates, strata = "treat", data = data_long_matching)
    print(tab1Unadj, test = FALSE, smd = TRUE)
    
    library(survey)
    cogSvy <- svydesign(ids = ~ HHIDPN, data = data_long_matching, weights = ~ mw)
    ## Weighted table with tableone
    tab1Mw <- svyCreateTableOne(vars = covariates, strata = "treat", data = cogSvy)
    print(tab1Mw, test = FALSE, smd = TRUE)
    
    ## Create SMD data frame
    dataPlot <- data.frame(variable = rownames(ExtractSmd(tab1Unadj)),
                           Unadjusted = ExtractSmd(tab1Unadj)[,"average"],
                           Weighted = ExtractSmd(tab1Mw)[,"average"])
    
    ## Reshape to long format
    library(reshape2)
    dataPlotMelt <- melt(data = dataPlot,
                         id.vars = "variable",
                         variable.name = "method",
                         value.name = "SMD")
    
    ## Variables names ordered by unadjusted SMD values
    varsOrderedBySmd <- rownames(dataPlot)[order(dataPlot[,"Unadjusted"])]
    varlabels<-c("Female","Years of education","Birth year","Father's education", "Mother's education", "Height","Age first married", "Age first smoked", "Age first worked","Never married", "Never smoked", "Never worked")
    
    ## Reorder factor levels
    dataPlotMelt$variable <- factor(dataPlotMelt$variable,
                                    #levels = varsOrderedBySmd)#,
                                    labels = varlabels)
    dataPlotMelt$method <- factor(dataPlotMelt$method,
                                  levels = c("Weighted","Unadjusted"),
                                  labels = c("Match weighted sample", "Unadjusted sample"))
    
    library(ggplot2)
    ggplot(data = dataPlotMelt, mapping = aes(x = variable, y = SMD, group = method)) +
      #geom_line() +
      geom_point(aes(colour = factor(method))) +
      scale_x_discrete(limits=rev)+
      geom_hline(yintercept = 0, size = 0.3) +
      geom_hline(yintercept = 0.1, size = 0.3) +
      #scale_x_discrete(labels = c("Female","Years of education","Birth year","Father's education", "Mother's education", "Height","Age first married", "Age first smoked", "Age first worked","Never married", "Never smoked", "Never worked"))+
      coord_flip() +
      theme_bw() + theme(legend.key = element_blank(),
                         axis.title.y=element_blank(),
                         legend.title=element_blank())
    
 #(B) Run regression analyses in match-weighted sample
    #Run match-weighted analyses
    table(data_long_matching$treat)
    data_long_matching = data_long_matching %>% 
      mutate(usbmxam = ifelse(treat == "dp_usbmxam", 1, 0),
             fbmxam = ifelse(treat == "dp_fbmxam", 1, 0))
    
    match.results = lmer(cog ~ AGE_R + factor(usbmxam) + factor(fbmxam) + 
                           AGE_R * factor(usbmxam) + AGE_R * factor(fbmxam) + 
                           factor(first_cog) + (AGE_R|HHIDPN), 
                         weights = mw, data = data_long_matching)
    
    match.unadj = summary(match.results)$coefficients
    (psm_coef<-match.unadj[c(3,4,6,7),1])
    psm_se<-match.unadj[c(3,4,6,7),2]
    (psm_lci<-psm_coef-(1.96*psm_se))
    (psm_uci<-psm_coef+(1.96*psm_se))
    
    #ShowRegTable(match.results,exp=FALSE)
    
