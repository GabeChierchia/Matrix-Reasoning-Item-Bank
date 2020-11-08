# Packages ---------------------------------------------------------------------------------------------------------
library(lmerTest)
library(emmeans)
library(doBy)
library(car)
library(Hmisc)
library(plyr)
library(xtable)
library(effects)
library(gridExtra)
library(moments)
library(data.table)
library(reshape2)
library(plotrix)
library(reshape2)
library(multicon)
library(merTools)
library(gtable)
library(lubridate)
library(ltm)



# Preprocessing---------------------------------------------------------------------------------------------------------
# Aggregate item-level data is available in the open science framework: 
# https://osf.io/g96f4/

# Code below reproduces results published in: 
#Chierchia, G., Fuhrmann, D., Knoll, L. J., Pi-Sunyer, B. P., Sakhardande, A.
#L., & Blakemore, S. J. (2019). The matrix reasoning item bank (MaRs-IB): novel,
#open-access abstract reasoning items for adolescents and adults. Royal Society
#open science, 6(10), 190232.

rr=read.csv("DIR", sep=",") # Load data 

# Relevel/relabel factors
rr$age4 = factor(rr$age4,levels(rr$age4)[c(4,2,3,1)])
rr$corrBin[rr$Correct==TRUE]=1
rr$corrBin[rr$Correct==FALSE]=0
rr$corrBin=as.factor(rr$corrBin)
rr$CorrectF=rr$corrBin
colnames(rr)[which(colnames(rr)=="Age")]="age"
# Standardize age
rr$age.s=scale(rr$age)

# Helmert coding factors
contrasts(rr$Relational.Reasoning..Distractor.Strategy) = contr.helmert(2)
contrasts(rr$puzzleSet) = contr.helmert(3)
contrasts(rr$Relational.Reasoning..Shape.Set) = contr.helmert(3)
contrasts(rr$age4) = contr.helmert(4)
contrasts(rr$Gender) = contr.helmert(2)

# Use below for figures
table <- summaryBy(Correct + Reaction.End ~ sub + School + age4 + age + Gender, rr, FUN=c(length,sum,mean,median))
tableCorr=summaryBy(Reaction.End ~ sub + age4 + age + Gender, droplevels(rr[rr$CorrectF=="1",]), FUN=c(median, IQR))

rr=droplevels(rr)

#--
# DESCRIPTIVES STATISTICS-------------------------------------------------------------------------------------------------------
#-
# DESCRIPTIVES OF TOTAL SAMPLE (REPORTED THROUGHOUT PAPER) -----
#-
# Not broken down by gender (used below) 
tableCorr2=summaryBy(Reaction.End ~ sub + age4, droplevels(rr[rr$CorrectF=="1",]), FUN=c(median, mean, sum))
colnames(tableCorr2)[c(3,4,5)]=c("median_RT_corr", "mean_RT_corr", "sum_RT_corr")
table2=merge(table, tableCorr2[,c(1,3,4,5)], all.x=T)
table2$ife_median=table2$median_RT_corr/table2$Correct.mean

# Overall stats
# Number of subs
nrow(table)

# Mean age
round(mean(table$age),2)
# Range age
round(range(table$age),2)
# SE age
round(std.error(table$age),2)
# SD age (reported on p. 6)
round(sd(table$age),2)

# reported on p. 11
# Mean correct 
round(mean(table$Correct.mean)*100, 2)
# SE correct
round(std.error(table$Correct.mean)*100, 2)
# Range3
round(range(table$Correct.mean)*100, 2)
# Grand median RTs (for items that were correctly answered)
round(median(tableCorr$Reaction.End.median), 2)
# IQR median rt correct
round(IQR(tableCorr$Reaction.End.median), 2)
# SE median rt correct
round(std.error(tableCorr$Reaction.End.median), 2)

# Trial number
trials = summaryBy(Trial.Number ~ sub + School + age4, rr, FUN=c(max))

# Mean, se, and range of trial numbers completed.
round(mean(trials$Trial.Number.max),2) # Reported on p. 11
round(sd(trials$Trial.Number.max)/sqrt(nrow(trials)),2)
range(trials$Trial.Number.max)

# Gender
xtabs(~Gender, table)

# What was the range of group size for adults? 
rr.a=rr[rr$age4=="Adults",]
#Ordering by subject, then by time stamp
rr.a=rr.a[order(rr.a$sub, as.Date(rr.a$Timestamp.Human)),]

# BY AGE AND GENDER (TABLE 1)-------------------------------------------------------
# inverse efficiency (IES)
IES.by.sub=table2
IES.by.age=data.frame(summaryBy(ife_median ~ age4, IES.by.sub, FUN=c(median, IQR)))
IES.by.age[,c(2,3)]=round(IES.by.age[,c(2,3)])

IES.by.age.and.gen=summaryBy(ife_median ~ age4 + Gender, table2, FUN=c(median, IQR))
IES.med.by.gen=data.frame(dcast(setDT(data.frame(IES.by.age.and.gen)), age4 ~ Gender, value.var="ife_median.median"))
IES.med.by.gen[,c(2,3)]=round(IES.med.by.gen[,c(2,3)])
IES.IQR.by.gen=data.frame(dcast(setDT(data.frame(IES.by.age.and.gen)), age4 ~ Gender, value.var="ife_median.IQR"))
IES.IQR.by.gen[,c(2,3)]=round(IES.IQR.by.gen[,c(2,3)])

IES.by.age$ife_median.median=paste0(IES.by.age$ife_median.median, " (m=", IES.med.by.gen$Male, " f=", IES.med.by.gen$Female, ")")
IES.by.age$ife_median.IQR=paste0(IES.by.age$ife_median.IQR, " (m=", IES.IQR.by.gen$Male, " f=", IES.IQR.by.gen$Female, ")")

rt.by.sub=summaryBy(Reaction.End ~ sub + age4 + age + Gender, droplevels(rr[rr$CorrectF=="1",]), FUN=c(median))
rt.by.age=data.frame(summaryBy(Reaction.End.median ~ age4, rt.by.sub, FUN=c(median, IQR)))
rt.by.age[,c(2,3)]=round(rt.by.age[,c(2,3)])

rt.by.age.and.gen=summaryBy(Reaction.End.median ~ age4 + Gender, rt.by.sub, FUN=c(median, IQR))
rt.med.by.gen=data.frame(dcast(setDT(data.frame(rt.by.age.and.gen)), age4 ~ Gender, value.var="Reaction.End.median.median"))
rt.med.by.gen[,c(2,3)]=round(rt.med.by.gen[,c(2,3)])
rt.IQR.by.gen=data.frame(dcast(setDT(data.frame(rt.by.age.and.gen)), age4 ~ Gender, value.var="Reaction.End.median.IQR"))
rt.IQR.by.gen[,c(2,3)]=round(rt.IQR.by.gen[,c(2,3)])

rt.by.age$Reaction.End.median.median=paste0(rt.by.age$Reaction.End.median.median, " (m=", rt.med.by.gen$Male, " f=", rt.med.by.gen$Female, ")")
rt.by.age$Reaction.End.median.IQR=paste0(rt.by.age$Reaction.End.median.IQR, " (m=", rt.IQR.by.gen$Male, " f=", rt.IQR.by.gen$Female, ")")

acc.by.sub <- summaryBy(Correct ~ sub + School + age4 + age + Gender, rr, FUN=c(length,mean))
acc.by.age=data.frame(summaryBy(Correct.mean + Correct.length ~ age4, acc.by.sub, FUN=c(mean, std.error, min, max)))
acc.by.age[,c(2:6)]=round(acc.by.age[,c(2:6)],2)
acc.by.age[,c(2,4,6,8)]=acc.by.age[,c(2,4,6,8)]*100

acc.by.age.and.gen=summaryBy(Correct.mean + Correct.length ~ age4 + Gender, acc.by.sub, FUN=c(mean, std.error, min, max))

acc.mean.mean.by.gen=dcast(setDT(acc.by.age.and.gen), age4 ~ Gender, value.var='Correct.mean.mean')
acc.mean.mean.by.gen[,c(2,3)]=100*round(acc.mean.mean.by.gen[,c(2,3)],2)
acc.by.age$Correct.mean.mean=paste0(acc.by.age$Correct.mean.mean, " (m=", acc.mean.mean.by.gen$Male, " f=", acc.mean.mean.by.gen$Female, ")")

acc.mean.se.by.gen=dcast(setDT(acc.by.age.and.gen), age4 ~ Gender, value.var='Correct.mean.std.error')
acc.mean.se.by.gen[,c(2,3)]=100*round(acc.mean.se.by.gen[,c(2,3)],2)
acc.by.age$Correct.mean.std.error=paste0(acc.by.age$Correct.mean.std.error, " (m=", acc.mean.se.by.gen$Male, " f=", acc.mean.se.by.gen$Female, ")")

acc.length.mean.by.gen=dcast(setDT(acc.by.age.and.gen), age4 ~ Gender, value.var='Correct.length.mean')
acc.length.mean.by.gen[,c(2,3)]=round(acc.length.mean.by.gen[,c(2,3)],2)
acc.by.age$Correct.length.mean=paste0(acc.by.age$Correct.length.mean, " (m=", acc.length.mean.by.gen$Male, " f=", acc.length.mean.by.gen$Female, ")")

acc.length.se.by.gen=dcast(setDT(acc.by.age.and.gen), age4 ~ Gender, value.var='Correct.length.std.error')
acc.length.se.by.gen[,c(2,3)]=round(acc.length.se.by.gen[,c(2,3)],2)
acc.by.age$Correct.length.std.error=paste0(acc.by.age$Correct.length.std.error, " (m=", acc.length.se.by.gen$Male, " f=", acc.length.se.by.gen$Female, ")")

acc.mean.max.by.gen=dcast(setDT(acc.by.age.and.gen), age4 ~ Gender, value.var='Correct.mean.max')
acc.mean.max.by.gen[,c(2,3)]=round(acc.mean.max.by.gen[,c(2,3)],2)*100
acc.by.age$Correct.mean.max=paste0(acc.by.age$Correct.mean.max, " (m=", acc.mean.max.by.gen$Male, " f=", acc.mean.max.by.gen$Female, ")")

acc.mean.min.by.gen=dcast(setDT(acc.by.age.and.gen), age4 ~ Gender, value.var='Correct.mean.min')
acc.mean.min.by.gen[,c(2,3)]=round(acc.mean.min.by.gen[,c(2,3)],2)*100
acc.by.age$Correct.mean.min=paste0(acc.by.age$Correct.mean.min, " (m=", acc.mean.min.by.gen$Male, " f=", acc.mean.min.by.gen$Female, ")")

acc.length.max.by.gen=dcast(setDT(acc.by.age.and.gen), age4 ~ Gender, value.var='Correct.length.max')
acc.by.age$Correct.length.max=paste0(acc.by.age$Correct.length.max, " (m=", acc.length.max.by.gen$Male, " f=", acc.length.max.by.gen$Female, ")")

acc.length.min.by.gen=dcast(setDT(acc.by.age.and.gen), age4 ~ Gender, value.var='Correct.length.min')
acc.by.age$Correct.length.min=paste0(acc.by.age$Correct.length.min, " (m=", acc.length.min.by.gen$Male, " f=", acc.length.min.by.gen$Female, ")")
summary(acc.by.age)
# Reordering columns: 
acc.by.age=acc.by.age[,c(1,2,4,6,8,3,5,7,9)]

descriptives=merge(rt.by.age, acc.by.age)
descriptives=merge(descriptives, IES.by.age)
# descriptives = descriptives[,c(1,4,8,2,3,6,10,12,13)] # There were errors here unfortunately. 
#colnames(descriptives)
descriptives = descriptives[,c('age4', "Correct.mean.mean", "Correct.mean.std.error", "Reaction.End.median.median", "Reaction.End.median.IQR", 
                               "Correct.length.mean", "Correct.length.std.error", "ife_median.median", "ife_median.IQR")] 
colnames(descriptives)=c("Age group", "Accuracy (mean)", "Accuracy (SE)", "RT (median)", "RT (IQR)", "Trials completed (mean)", 
                         "Trials completed (SE)", "Inverse efficiency (median)", "Inverse efficiency (IQR)")
descriptives=t(descriptives)

descriptives=descriptives[,c(4,2,3,1)]

write.csv(descriptives,paste0("DIR", Sys.Date(),".txt"))
write.csv(descriptives, "DIR")


#-
# BY PUZZLE SET (SUPPLEMENTARY TABLE 1)  ------------------------------------------------------
#-

rt.by.sub=summaryBy(Reaction.End ~ sub + puzzleSet, droplevels(rr[rr$CorrectF=="1" & !is.na(rr$puzzleSet),]), FUN=c(median))
rt.by.age=data.frame(summaryBy(Reaction.End.median ~ puzzleSet, rt.by.sub, FUN=c(median, IQR)))
rt.by.age[,c(2,3)]=round(rt.by.age[,c(2,3)])

acc.by.sub <- summaryBy(Correct ~ sub + puzzleSet, rr[!is.na(rr$puzzleSet),], FUN=c(length,mean))
acc.by.age=data.frame(summaryBy(Correct.mean + Correct.length ~ puzzleSet, acc.by.sub, FUN=c(mean, std.error, min, max)))
acc.by.age[,c(2:6)]=round(acc.by.age[,c(2:6)],2)
acc.by.age[,c(2,4,6,8)]=acc.by.age[,c(2,4,6,8)]*100

acc.by.puzzleSet=acc.by.age[,c(1,2,4,6,8,3,5,7,9)]

descriptives=merge(rt.by.age, acc.by.puzzleSet)
descriptives=t(descriptives)
#descriptives=descriptives[,c(4,2,3,1)]
write.csv(descriptives, "DIR")

#-
# Desctractor strategy (SUPPLEMENTARY TABLE 2)  ------------------------------------------------------
#-

rt.by.sub=summaryBy(Reaction.End ~ sub + Relational.Reasoning..Distractor.Strategy, droplevels(rr[rr$CorrectF=="1" & !is.na(rr$puzzleSet),]), FUN=c(median))
rt.by.age=data.frame(summaryBy(Reaction.End.median ~Relational.Reasoning..Distractor.Strategy, rt.by.sub, FUN=c(median, IQR)))
rt.by.age[,c(2,3)]=round(rt.by.age[,c(2,3)])

acc.by.sub <- summaryBy(Correct ~ sub + Relational.Reasoning..Distractor.Strategy, rr[!is.na(rr$puzzleSet),], FUN=c(length,mean))
acc.by.age=data.frame(summaryBy(Correct.mean + Correct.length ~Relational.Reasoning..Distractor.Strategy, acc.by.sub, FUN=c(mean, std.error, min, max)))
acc.by.age[,c(2:6)]=round(acc.by.age[,c(2:6)],2)
acc.by.age[,c(2,4,6,8)]=acc.by.age[,c(2,4,6,8)]*100

acc.by.distractor=acc.by.age[,c(1,2,4,6,8,3,5,7,9)]

descriptives=merge(rt.by.age, acc.by.distractor)
descriptives=t(descriptives)
#descriptives=descriptives[,c(4,2,3,1)]
write.csv(descriptives, "DIR")

#-
# Shape set (SUPPLEMENTARY TABLE 3) ------------------------------------------------------
#-

rt.by.sub=summaryBy(Reaction.End ~ sub + Relational.Reasoning..Shape.Set, droplevels(rr[rr$CorrectF=="1" & !is.na(rr$puzzleSet),]), FUN=c(median))
rt.by.age=data.frame(summaryBy(Reaction.End.median ~ Relational.Reasoning..Shape.Set, rt.by.sub, FUN=c(median, IQR)))
rt.by.age[,c(2,3)]=round(rt.by.age[,c(2,3)])

acc.by.sub <- summaryBy(Correct ~ sub + Relational.Reasoning..Shape.Set, rr[!is.na(rr$puzzleSet),], FUN=c(length,mean))
acc.by.age=data.frame(summaryBy(Correct.mean + Correct.length ~Relational.Reasoning..Shape.Set, acc.by.sub, FUN=c(mean, std.error, min, max)))
acc.by.age[,c(2:6)]=round(acc.by.age[,c(2:6)],2)
acc.by.age[,c(2,4,6,8)]=acc.by.age[,c(2,4,6,8)]*100

acc.by.shape.set=acc.by.age[,c(1,2,4,6,8,3,5,7,9)]

descriptives=merge(rt.by.age, acc.by.shape.set)
descriptives=t(descriptives)
write.csv(descriptives, "DIR")

#-
# Practice trials info--------------------------------------------------------------------------------------------------------
#-
practice=read.csv("DIR", sep=',')

colnames(practice)[1]="sub"
sumfun <- function(x, ...){
  c(sum=sum(x,na.rm=T),length=length(x))
}
table.practice <- summaryBy(Correct ~ sub + School + age4, practice, FUN=sumfun)
table.repeaters=droplevels(table.practice[table.practice$Correct.length >3,])
# How many ppts needed 3 or more? 
nrow(table.repeaters)
# Of those who did, average attempts? 
round(mean(table.repeaters$Correct.length),2)
round(std.error(table.repeaters$Correct.length),2)
# Never more than? 
max(table.repeaters$Correct.length)

#-
# Skew, kurtosis, floor and ceiling performance----------------------------------------------------------------
#-

# Mean accuracy
round(mean(table$Correct.mean)*100, 2)
round(std.error(table$Correct.mean)*100, 2)
round(range(table$Correct.mean)*100, 2)
round(skewness(table$Correct.mean),2)
round(kurtosis(table$Correct.mean),2)
floor = length(table$Correct.mean[table$Correct.mean <= 0.25])
floor
ceiling = length(table$Correct.mean[table$Correct.mean > 0.9999])
ceiling
round((1 - (floor+ceiling)/length(table$Correct.mean))*100,2)

# split half reliability (Spearman - Brown prophecy formula)
temp = subset(rr, select=c("sub","Trial.Number","corrBin")) 
TrialData = dcast(temp, sub ~ Trial.Number, value.var="corrBin")
TrialData = lapply(TrialData[-1], as.character) # delete sub column, convert to char
TrialData = data.frame(lapply(TrialData, as.numeric)) # make numeric
round(splithalf.r(TrialData, sims = 1000, graph = T, seed = 2),2)

# For each age group separately: 
by.age=data.frame(age=levels(rr$age4), mean.split.half.r=NA)
for (i in 1:nrow(by.age)){
temp = subset(droplevels(rr[rr$age4==by.age$age[i],]), select=c("sub","Trial.Number","corrBin")) 
TrialData = dcast(temp, sub ~ Trial.Number, value.var="corrBin")
TrialData = lapply(TrialData[-1], as.character) # delete sub column, convert to char
TrialData = data.frame(lapply(TrialData, as.numeric)) # make numeric
test=splithalf.r(TrialData, sims = 1000, graph = T, seed = 2)
by.age$mean.split.half.r[i]=round(test[3],2)
}
by.age
# Test-retest
trt=read.csv("DIR")
trt=droplevels(trt[trt$group==levels(trt$group)[2],])

# Calculating days elapsed between t1 and t2: I have two dates for some reason,
# I checked the second however and it seems to be invariant for each
# participant, therefore it is probably unrelated to the questions of interest
# here.
trt$diff.t2.t1.x = as.numeric(as.Date(as.character(trt$date.x.2), format="%d/%m/%Y")-
                              as.Date(as.character(trt$date.x.1), format="%d/%m/%Y"))

# trt$diff.t2.t1.y = as.numeric(as.Date(as.character(trt$date.y.2), format="%Y-%m-%d")-
#                                 as.Date(as.character(trt$date.y.1), format="%Y-%m-%d"))

# How many ppts? 
length(unique(na.omit(trt$id)))

# Stats on days elapsed
range(trt$diff.t2.t1.x, na.rm=T)
round(std.error(trt$diff.t2.t1.x, na.rm=T))
round(mean(trt$diff.t2.t1.x, na.rm=T))

cor.test(trt[,"t1"], trt[,"t2"]) # Correlation reported. 
contrasts(trt$age4) = contr.helmert(4)
t.ret=lmer(t1 ~ age4 * t2 + (1|school), trt)
Anova(t.ret, type=3) # Absence of interaction with age reported. 

# DIGIT-SPAN TO RR ()
data=read.table("DIR", header=T)
dat_DS <- subset(data, 
                 Section=="digit_span_trials" & Screen.Type=="trial" & Reaction.End>250 & Correct=="TRUE" & Procedure.Stage=="T1",
                 select = c(Participant.ID, Group.x,Screen.Type,Trial.Number,
                            Procedure.Stage, Correct,Section,Difficulty.Raw))
dat_corr <- summaryBy(Difficulty.Raw~Participant.ID+Group.x,data=dat_DS, FUN=max)
colnames(dat_corr)[1]="sub"

table <- summaryBy(Correct + Reaction.End ~ sub + School + age4 + age + Gender, rr, FUN=c(length,sum,mean,median))
ds.rr=merge(table, dat_corr, all.x=T)

# Did the measures correlate? 
with(ds.rr, cor.test(Correct.mean, Difficulty.Raw.max))
with(ds.rr, cor.test(Correct.mean, Difficulty.Raw.max,method="spearman"))
# Yes, moderately. 

# Were there interactions with age?
contrasts(table$age4) = contr.helmert(4)
ds.mod=lmer(Correct.mean ~ age4 * Difficulty.Raw.max + (1|School), ds.rr)
Anova(ds.mod, type=3) 
# No. 

####################################################################################################################################################
# EFFECTS OF DISTRACTOR STRATEGY, SHAPE SET, PUZZLE SET & DIFFICULTY (Mixed models)
####################################################################################################################################################

############
# ACCURACY
############
model_corr_distractor <- glmer(CorrectF~ 
                                 Relational.Reasoning..Distractor.Strategy + (1|School/sub)+ (1|barcode), 
                               data=rr, 
                               family = binomial,
                               control = glmerControl(optimizer="bobyqa",
                                                      optCtrl=list(maxfun=100000)))
Anova(model_corr_distractor, type="III") # Reported

model_corr_puzzle <- glmer(CorrectF~ 
                             puzzleSet + (1|School/sub) + (1|barcode), 
                           data=rr, 
                           family = binomial,
                           control = glmerControl(optimizer="bobyqa",
                                                  optCtrl=list(maxfun=100000)))
Anova(model_corr_puzzle, type="III") # Reported

model_corr_shape <- glmer(CorrectF~ 
                            Relational.Reasoning..Shape.Set + (1|School/sub)+ (1|barcode), 
                          data=rr, 
                          family = binomial,
                          control = glmerControl(optimizer="bobyqa",
                                                 optCtrl=list(maxfun=100000)))
Anova(model_corr_shape, type="III") # Reported

model_corr_diff <- glmer(CorrectF~ 
                            diff + (1|School/sub), 
                          data=rr, 
                          family = binomial,
                          control = glmerControl(optimizer="bobyqa",
                                                 optCtrl=list(maxfun=100000)))
Anova(model_corr_diff, type="III") # Reported
summary(model_corr_diff)

table.diff <- summaryBy(Correct + Reaction.End ~ sub + diff, rr, FUN=c(length,sum,mean,median))
cor.test(table.diff$Correct.mean, table.diff$diff)
cor.test(table.diff$Reaction.End.median, table.diff$diff)


############
# RTs
############
model_RT_puzzle  <- glmer(Reaction.End~ 
                        puzzleSet +(1|School/sub)+ (1|barcode), 
                      data=rr[rr$CorrectF==1 & !is.na(rr$Relational.Reasoning..Shape.Set),],
                      family=Gamma(link="log"))
Anova(model_RT_puzzle, type="III")# Reported

model_RT_shape  <- glmer(Reaction.End~ 
                          Relational.Reasoning..Shape.Set +(1|School/sub)+ (1|barcode), 
                        data=rr[rr$CorrectF==1 & !is.na(rr$Relational.Reasoning..Shape.Set),],
                        family=Gamma(link="log"))
Anova(model_RT_shape, type="III") # Reported
emmeans(model_RT_shape, pairwise ~ Relational.Reasoning..Shape.Set, adjust="Bonferroni", type="response")

model_RT_distractor  <- glmer(Reaction.End~ 
                        Relational.Reasoning..Distractor.Strategy +(1|School/sub)+ (1|barcode), 
                      data=rr[rr$CorrectF==1& !is.na(rr$Relational.Reasoning..Shape.Set),],
                      family=Gamma(link="log"))
Anova(model_RT_distractor, type="III") # Reported

model_RT_diff  <- glmer(Reaction.End~ 
                            diff +(1|School/sub), 
                          data=rr[rr$CorrectF==1 & !is.na(rr$Relational.Reasoning..Shape.Set),],
                          family=Gamma(link="log"))
Anova(model_RT_diff, type="III")# Reported
summary(model_RT_diff)
#-
# EFFECTS OF AGE AND GENDER ON ACCURACY-----------------------------------------------------------------
#-


#-
# AGE DISCRETE
#-
# Models
age_model_corr <- glmer(CorrectF~ 
                        age4 + (1|School/sub)+ (1|barcode), 
                        data=rr, 
                        family = binomial,
                        control = glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=100000)))
#qqnorm(residuals(age_model_corr))# Mmm

Anova(age_model_corr, type="III") # Reported
emmeans(age_model_corr, pairwise ~ age4,adjust="Bonferroni") # Reported
plot(allEffects(age_model_corr))

# Tables
cons=as.data.frame(emmeans(age_model_corr, pairwise ~ age4, adjust="Bonferroni")$contrasts)
cons$p.uncorrected=round(as.data.frame(emmeans(age_model_corr, pairwise ~ age4, adjust="none")$contrasts)$p.value,3)
cons$estimate=round(cons$estimate,2)
cons$SE=round(cons$SE,2)
cons$p.value=round(cons$p.value,3)
cons=cons[,-which(colnames(cons)%in% c("z.ratio", "df"))]
write.csv(cons, "/Users/gabrielechierchia/Desktop/Blakemore/relational reasoning/tables/contrasts in accuracy by age group.txt", row.names=F)

# Plots
# Method for extracting estimates partially documented here:
# https://danieljhocking.wordpress.com/2012/07/25/plotting-95-confidence-bands-in-r/
m.newdat <- expand.grid(
  age4=levels(rr$age4), 
  CorrectF = 0
) 
contrasts(m.newdat$age4) = contr.helmert(4)
m.mm <- model.matrix(terms(age_model_corr),m.newdat) 
m.newdat$CorrectF <- m.mm %*% fixef(age_model_corr) 
m.pvar1 <- diag(m.mm %*% tcrossprod(vcov(age_model_corr),m.mm)) 
m.newdat <- data.frame(
  m.newdat
  , m.plo = plogis(m.newdat$CorrectF-2*sqrt(m.pvar1)) 
  , m.phi = plogis(m.newdat$CorrectF+2*sqrt(m.pvar1))
)

m.newdat$CorrectF=plogis(m.newdat$CorrectF) # back-transforming predicted probabilities

# The two lines below show that the manually pulled means are identical to those returned by emmeans: 
# emmeans.prob=as.data.frame(emmeans(age_model_corr, pairwise ~ age4,adjust="Bonferroni", type="response")$emmeans)$prob
# man.pull.prob=as.numeric(m.newdat$CorrectF)
# cbind(man.pull.prob, emmeans.prob)

# The confidence intervals are only slightly off (max by 1 2nd decimal digits,
# possibly due to different degree of freedom estimation?)

ggplot(table, aes(x=age4, y=Correct.mean, fill=age4)) +
  geom_violin(alpha=0.8)+theme_bw()+#geom_boxplot(width=0.08, fill='white')+
  geom_point(position=position_jitter(width=0.1, height=0), size=1, alpha=0.2)+
  ylab("Accuracy")+xlab("Age group")+
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=22,face="bold"), 
        panel.grid.minor = element_blank())+
  guides(fill=FALSE) +
  geom_point(data=m.newdat, aes(age4, CorrectF), shape=15, size=3)+
  geom_errorbar(data=m.newdat, aes(age4, CorrectF, ymin = m.plo, ymax = m.phi),width = 0.05,size  = 0.5)+
  coord_cartesian(ylim = c(0, 1)) 

ggsave('DIR', width=12, height=8)

#-
# AGE CONTINUOUS
#-
# Generate orthogonolized polynomial components 
rr$age.1=poly(rr$age,3)[,1]
rr$age.2=poly(rr$age,3)[,2]
rr$age.3=poly(rr$age,3)[,3]

# Models
age_model_corr_cont <- glmer(CorrectF~ 
                               poly(age,3) + (1|School/sub) + (1|barcode), 
                             data=rr, # Replace with rr for full data set, rr.red.3 for reduced dataset (i.e., no variance in completion rates)
                             family = binomial,
                             control = glmerControl(optimizer="bobyqa",
                                                    optCtrl=list(maxfun=100000)))

Anova(age_model_corr_cont, type="III") # Reported
anova(age_model_corr_cont, type="III") # Reported
summary(age_model_corr_cont)

# The model below is only to get separate chisq for each age poly (which are not
# given above when using poly(age,3) directly in model specification.
age_model_corr_cont2 <- glmer(CorrectF~
                                age.1 + age.2 + age.3 + (1|School/sub) + (1|barcode),
                              data=rr,
                              family = binomial,
                              control = glmerControl(optimizer="bobyqa",
                                                     optCtrl=list(maxfun=100000)))
 Anova(age_model_corr_cont2, type="III") # Reported

# Plots
m.newdat <- expand.grid(
  age=unique(rr$age), 
  CorrectF= 0
) 
m.mm <- model.matrix(terms(age_model_corr_cont),m.newdat) 
m.newdat$CorrectF <- m.mm %*% fixef(age_model_corr_cont) 
m.pvar1 <- diag(m.mm %*% tcrossprod(vcov(age_model_corr_cont),m.mm)) 
m.newdat <- data.frame(
  m.newdat
  , m.plo = plogis(m.newdat$CorrectF-2*sqrt(m.pvar1))
  , m.phi = plogis(m.newdat$CorrectF+2*sqrt(m.pvar1))
)

m.newdat$CorrectF=plogis(m.newdat$CorrectF)

ggplot(table, aes(age, Correct.mean))+geom_point(alpha=0.3, size=0.8)+geom_line(data=m.newdat,aes(age, CorrectF), size=1) + 
  geom_ribbon(data=m.newdat, aes(age,CorrectF, ymin=m.plo, ymax=m.phi),alpha=0.2)+theme_bw()+
  xlab("Age")+ylab("")+
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=22,face="bold"),
        panel.grid.minor = element_blank())+
  coord_cartesian(ylim = c(0, 1)) 

ggsave('DIR',width=10, height=8)

#-
# Components separately
#-
# Plots
m.newdat <- expand.grid(
  age=unique(rr$age), 
  CorrectF= 0
) 

m.mm <- model.matrix(terms(age_model_corr_cont),m.newdat) 
m.newdat$CorrectF <- m.mm %*% fixef(age_model_corr_cont) 
m.pvar1 <- diag(m.mm %*% tcrossprod(vcov(age_model_corr_cont),m.mm)) 
m.newdat <- data.frame(
  m.newdat
  , m.plo = plogis(m.newdat$CorrectF-2*sqrt(m.pvar1))
  , m.phi = plogis(m.newdat$CorrectF+2*sqrt(m.pvar1))
)

m.newdat$CorrectF=plogis(m.newdat$CorrectF)
m.newdat$CorrectF.1 <- plogis(m.mm[,c(1,2)] %*% fixef(age_model_corr_cont)[c(1,2)])
m.newdat$CorrectF.2  <- plogis(m.mm[,c(1,3)] %*% fixef(age_model_corr_cont)[c(1,3)])
m.newdat$CorrectF.3  <- plogis(m.mm[,c(1,4)] %*% fixef(age_model_corr_cont)[c(1,4)])

ggplot(table, aes(age, Correct.mean))+geom_point(alpha=0.3, size=0.8)+theme_bw()+
  xlab("Age")+ylab("")+geom_line(data=m.newdat,aes(age, CorrectF), size=1) + 
  geom_ribbon(data=m.newdat, aes(age,CorrectF, ymin=m.plo, ymax=m.phi),alpha=0.2)+theme_bw()+
  geom_line(data=m.newdat,aes(x=age, y=CorrectF.1, col='1linear'), size=0.8, linetype=2)+
  geom_line(data=m.newdat,aes(x=age, y=CorrectF.2, col='2quadratic'), size=0.8, linetype=2)+
  #geom_line(data=m.newdat,aes(x=age, y=CorrectF.3, col='3cubic'), size=0.8, linetype=2)+theme_bw()+
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=22,face="bold"),
        panel.grid.minor = element_blank(),
        legend.position="none")+coord_cartesian(ylim = c(0, 1))

ggsave('DIR',width=12, height=8)
ggsave('DIR',width=12, height=8) # Added during proofs


#-
# GENDER
#-
gender_model_corr <- glmer(CorrectF~ 
                             age4 * Gender + (1|School/sub)+ (1|barcode), 
                           data=rr, 
                           family = binomial,
                           control = glmerControl(optimizer="bobyqa",
                                                  optCtrl=list(maxfun=100000)))
Anova(gender_model_corr, type="III") # Reported
summary(gender_model_corr)
emmeans(gender_model_corr, pairwise ~ Gender,adjust="Bonferroni") # Reported


# Plot
gen.acc=data.frame(emmeans(gender_model_corr, pairwise ~ Gender,adjust="Bonferroni", type="response")$emmeans)
ggplot(table, aes(x=Gender, y=Correct.mean, fill=Gender)) +
  geom_violin(alpha=0.8)+theme_bw()+#geom_boxplot(width=0.08, fill='white')+
  geom_point(position=position_jitter(width=0.1, height=0), size=1, alpha=0.2)+ylab("Accuracy")+xlab("")+
  theme(axis.text=element_text(size=28),axis.title=element_text(size=28,face="bold"))+
  scale_fill_discrete(name = "Gender")+coord_cartesian(ylim = c(0, 1))+guides(fill=FALSE) +
  geom_point(data=gen.acc, aes(Gender, prob), shape=15, size=3)+
  geom_errorbar(data=gen.acc, aes(Gender, prob, ymin = asymp.LCL, ymax = asymp.UCL),width = 0.05,size  = 0.5)
ggsave('DIR',width=10, height=8)

#-
# EFFECTS OF AGE AND GENDER ON RTs----------------------------------------------------------------------
#-

#-
# AGE DISCRETE
#-
age_model_RT <- glmer(Reaction.End~ 
                        age4 + (1|School/sub)+ (1|barcode), 
                       data=rr[rr$CorrectF=="1",],
                     family=Gamma(link="log"))
Anova(age_model_RT, type="III")
emmeans(age_model_RT, pairwise ~ age4,adjust="none", type="response") 
emmeans(age_model_RT, pairwise ~ age4,adjust="none", type="response") 
emmeans(age_model_RT, pairwise ~ age4,adjust="Bonferroni", type="response") # Reported
emmeans(age_model_RT, pairwise ~ age4,adjust="Bonferroni")

age_model_RT.log <- lmer(log(Reaction.End)~ 
                        age4 + (1|School/sub)+ (1|barcode), 
                      data=rr[rr$CorrectF=="1",])
Anova(age_model_RT.log, type="III")

# Tables
cons=as.data.frame(emmeans(age_model_RT, pairwise ~ age4, adjust="Bonferroni")$contrasts)
cons$p.uncorrected=round(as.data.frame(emmeans(age_model_RT, pairwise ~ age4, adjust="none")$contrasts)$p.value,3)
cons$estimate=round(cons$estimate,2)
cons$SE=round(cons$SE,2)
cons$p.value=round(cons$p.value,3)
cons=cons[,-which(colnames(cons)%in% c("z.ratio", "df"))]
write.csv(cons, "DIR", row.names=F)

# Plotting
m.newdat <- expand.grid(
  age4=levels(rr$age4), 
  Reaction.End = 0
) 
contrasts(m.newdat$age4) = contr.helmert(4)
m.mm <- model.matrix(terms(age_model_RT),m.newdat) 
m.newdat$Reaction.End <- m.mm %*% fixef(age_model_RT) 
m.pvar1 <- diag(m.mm %*% tcrossprod(vcov(age_model_RT),m.mm)) 
m.newdat <- data.frame(
  m.newdat
  , m.plo = exp(m.newdat$Reaction.End-2*sqrt(m.pvar1))
  , m.phi = exp(m.newdat$Reaction.End+2*sqrt(m.pvar1))
)

m.newdat$Reaction.End=exp(m.newdat$Reaction.End) # back-transforming predicted probabilities

ggplot(tableCorr, aes(x=age4, y=Reaction.End.median, fill=age4)) +
  geom_violin(alpha=0.8)+theme_bw()+#geom_boxplot(width=0.08, fill='white')+
  geom_point(position=position_jitter(width=0.1, height=0), size=1, alpha=0.2)+ylab("Response times")+xlab("Age Group")+
  theme(axis.text=element_text(size=22),axis.title=element_text(size=22,face="bold"))+
  guides(fill=FALSE) +
  geom_point(data=m.newdat, aes(age4, Reaction.End), shape=15, size=3)+
  geom_errorbar(data=m.newdat, aes(age4, Reaction.End, ymin = m.plo, ymax = m.phi),width = 0.05,size  = 0.5)
ggsave('DIR',width=10, height=8)

#-
# AGE CONTINUOUS
#-
# Modelling
age_model_rt_cont <- glmer(Reaction.End ~ 
                             poly(age,3) + (1|School/sub)+ (1|barcode), 
                           data=rr[rr$CorrectF=="1",],
                           family=Gamma(link="log"), control = glmerControl(optimizer="bobyqa",
                                                                            optCtrl=list(maxfun=100000)))

# 
 age_model_rt_cont.2 <- glmer(Reaction.End ~ 
                                age.1 + age.2 + (1|School/sub)+ (1|barcode), 
                              data=rr.red3[rr.red3$CorrectF=="1",],
                              family=Gamma(link="log"),
                              control = glmerControl(optimizer="bobyqa",
                                                     optCtrl=list(maxfun=100000)))
Anova(age_model_rt_cont, type="III") 
summary(age_model_rt_cont) # curious, the summary shows a significant linear trend. 

age_model_rt_cont.log <- lmer(log(Reaction.End) ~ 
                             poly(age,3) + (1|School/sub)+ (1|barcode), 
                           data=rr[rr$CorrectF=="1",])
Anova(age_model_rt_cont.log,type=3)
summary(age_model_rt_cont.log)

# Dropping the cubic component allows for convergence
age_model_rt_cont.log <- lmer(log(Reaction.End) ~ 
                                age.1 +age.2+ (1|School/sub)+ (1|barcode), 
                              data=rr.red3[rr.red3$CorrectF=="1",])
Anova(age_model_rt_cont.log,type=3)
summary(age_model_rt_cont.log)
AIC(age_model_rt_cont.log)
AIC(age_model_rt_cont.2)

# Plotting
m.newdat <- expand.grid(
  age=unique(rr$age), 
  Reaction.End= 0
) 

m.mm <- model.matrix(terms(age_model_rt_cont),m.newdat) 
m.newdat$Reaction.End <- m.mm %*% fixef(age_model_rt_cont) 
m.pvar1 <- diag(m.mm %*% tcrossprod(vcov(age_model_rt_cont),m.mm)) 
m.newdat <- data.frame(
  m.newdat
  , m.plo = exp(m.newdat$Reaction.End-2*sqrt(m.pvar1))
  , m.phi = exp(m.newdat$Reaction.End+2*sqrt(m.pvar1))
)

m.newdat$Reaction.End=exp(m.newdat$Reaction.End)

ggplot(tableCorr, aes(age,Reaction.End.median))+ geom_point(alpha=0.3, size=0.8) + 
  geom_line(data=m.newdat,aes(x=age, y=Reaction.End), size=1)+
  geom_ribbon(data=m.newdat, aes(x=age,y=Reaction.End, ymin=m.plo, ymax=m.phi), alpha=0.2)+theme_bw()+
  xlab("Age")+ylab("")+
  theme(axis.text=element_text(size=22),axis.title=element_text(size=22,face="bold"))
ggsave('DIR',width=10, height=8)


#-
# Components separately
#-
m.newdat <- expand.grid(
  age=unique(rr$age), 
  Reaction.End= 0
) 

m.mm <- model.matrix(terms(age_model_rt_cont),m.newdat) 
m.newdat$Reaction.End <- m.mm %*% fixef(age_model_rt_cont) 
m.pvar1 <- diag(m.mm %*% tcrossprod(vcov(age_model_rt_cont),m.mm)) 
m.newdat <- data.frame(
  m.newdat
  , m.plo = exp(m.newdat$Reaction.End-2*sqrt(m.pvar1))
  , m.phi = exp(m.newdat$Reaction.End+2*sqrt(m.pvar1))
)

m.newdat$Reaction.End=exp(m.newdat$Reaction.End)
m.newdat$Reaction.End.1 <- exp(m.mm[,c(1,2)] %*% fixef(age_model_rt_cont)[c(1,2)])
m.newdat$Reaction.End.2  <- exp(m.mm[,c(1,3)] %*% fixef(age_model_rt_cont)[c(1,3)])
m.newdat$Reaction.End.3  <- exp(m.mm[,c(1,4)] %*% fixef(age_model_rt_cont)[c(1,4)])

# This just manually picks up ggplot default colors. I need them because I want 
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

ggplot(tableCorr, aes(age,Reaction.End.median))+ geom_point(alpha=0.3, size=0.8) +
  geom_line(data=m.newdat,aes(x=age, y=Reaction.End), size=1)+
  geom_ribbon(data=m.newdat, aes(x=age,y=Reaction.End, ymin=m.plo, ymax=m.phi), alpha=0.2)+
  geom_line(data=m.newdat,aes(x=age, y=Reaction.End.1, colour='red'), size=0.8, linetype=2)+theme_bw()+
  xlab("Age")+ylab("")+
  theme(axis.text=element_text(size=22),axis.title=element_text(size=22,face="bold"),legend.position="none")
ggsave('DIR',width=10, height=8)


#-
# GENDER
#-

# Model variants with log(rt) as dependent variable (with gamma and
# log-link and without log-link). They did not converge. 
gender_model_RT <- lmer(log(Reaction.End)~ 
                       Gender* age4 + (1|School/sub)+ (1|barcode), 
                     data=rr[rr$CorrectF=="1",])
Anova(gender_model_RT, type="III")

emmeans(gender_model_RT, pairwise ~ Gender,adjust="Bonferroni", type="response") 
emmeans(gender_model_RT, pairwise ~ Gender,adjust="Bonferroni") 

gen.rt=data.frame(emmeans(gender_model_RT, pairwise ~ Gender,adjust="Bonferroni", type="response")$emmeans)
# Gender plot
ggplot(tableCorr, aes(x=Gender, y=Reaction.End.median, fill=Gender)) +
  geom_violin(alpha=0.8)+theme_bw()+#geom_boxplot(width=0.08, fill='white')+
  geom_point(position=position_jitter(width=0.1, height=0), size=1, alpha=0.2)+ylab("Response times")+xlab("")+
  theme(axis.text=element_text(size=28),axis.title=element_text(size=28,face="bold"))+
  guides(fill=FALSE) +
  geom_point(data=gen.rt, aes(Gender, response), shape=15, size=3)+
  geom_errorbar(data=gen.rt, aes(Gender, response, ymin = asymp.LCL, ymax = asymp.UCL),width = 0.05,size  = 0.5)
ggsave('DIR',width=10, height=8)

#-
# EFFECTS OF AGE AND GENDER ON SPEED/ACCURACY----------------------------------------------------------------------
#-
# Choose one: 
tableCorr2=summaryBy(Reaction.End ~ sub + age4, droplevels(rr[rr$CorrectF=="1",]), FUN=c(median, mean, sum))

colnames(tableCorr2)[c(3,4,5)]=c("median_RT_corr", "mean_RT_corr", "sum_RT_corr")

table2=merge(table, tableCorr2[,c(1,3,4,5)], all.x=T)
table2$age.1=poly(table2$age,3)[,1]
table2$age.2=poly(table2$age,3)[,2]
table2$age.3=poly(table2$age,3)[,3]

contrasts(table2$age4) = contr.helmert(4)
contrasts(table2$Gender) = contr.helmert(2)

# Bruyer & BRYSBAERT 2011
# Ife=mean (or median) RT of correct responses divided by proportion of correct responses
table2$ife_median=table2$median_RT_corr/table2$Correct.mean

#-
# AGE DISCRETE
#-
# Modeling
ife_age=lmer(ife_median ~ age4 + (1|School), table2)
Anova(ife_age, type="III") 
emmeans(ife_age, pairwise ~ age4, adjust="Bonferroni") 
emmeans(ife_age, pairwise ~ age4, adjust="none")

# Tables
cons=as.data.frame(emmeans(ife_age, pairwise ~ age4, adjust="Bonferroni")$contrasts)
cons$p.uncorrected=round(as.data.frame(emmeans(ife_age, pairwise ~ age4, adjust="none")$contrasts)$p.value,3)
cons$estimate=round(cons$estimate,2)
cons$SE=round(cons$SE,2)
cons$p.value=round(cons$p.value,3)
cons=cons[,-which(colnames(cons)%in% c("t.ratio", "df"))]
write.csv(cons, "DIR", row.names=F)


# Plots
m.newdat <- expand.grid(
  age4=levels(table2$age4), 
  ife_median = 0
) 
#contrasts(m.newdat$age4) = contr.helmert(4)
m.mm <- model.matrix(terms(ife_age),m.newdat) 
m.newdat$ife_median <- m.mm %*% fixef(ife_age) 
m.pvar1 <- diag(m.mm %*% tcrossprod(vcov(ife_age),m.mm)) 
m.newdat <- data.frame(
  m.newdat
  , m.plo = m.newdat$ife_median-2*sqrt(m.pvar1)
  , m.phi = m.newdat$ife_median+2*sqrt(m.pvar1)
  #, m.tlo = m.newdat$respC-2*sqrt(m.tvar1)
  #, m.thi = m.newdat$respC+2*sqrt(m.tvar1)
)
dev.off()
ggplot(table2, aes(x=age4, y=ife_median, fill=age4)) +
  geom_violin(alpha=0.8)+theme_bw()+#geom_boxplot(width=0.08, fill='white')+
  geom_point(position=position_jitter(width=0.1, height=0), size=1, alpha=0.2)+ylab("Inverse efficiency")+xlab("Age Group")+
  theme(axis.text=element_text(size=22),axis.title=element_text(size=22,face="bold"))+
  scale_fill_discrete(name = "Age group")+guides(fill=FALSE) +
  geom_point(data=m.newdat, aes(age4, ife_median), shape=15, size=3)+
  geom_errorbar(data=m.newdat, aes(age4, ife_median, ymin = m.plo, ymax = m.phi),width = 0.05,size  = 0.5)
ggsave('DIR',width=10, height=8)

#-
# AGE CONTINUOUS
#-
ife_age_cont = lmer(ife_median ~ poly(age,3)  + (1|School), table2)
Anova(ife_age_cont, type="III") 
summary(ife_age_cont)
summary(ife_age_cont,ddf="Kenward-Roger")

m.newdat <- expand.grid(
  age=unique(rr$age), 
  ife_median= 0
) 

m.mm <- model.matrix(terms(ife_age_cont),m.newdat) 
m.newdat$ife_median <- m.mm %*% fixef(ife_age_cont) 
m.pvar1 <- diag(m.mm %*% tcrossprod(vcov(ife_age_cont),m.mm)) 
m.newdat <- data.frame(
  m.newdat
  , m.plo = m.newdat$ife_median-2*sqrt(m.pvar1)
  , m.phi = m.newdat$ife_median+2*sqrt(m.pvar1)
)

ggplot(table2, aes(age,ife_median))+ geom_point(alpha=0.3, size=0.8) + 
  geom_line(data=m.newdat,aes(x=age, y=ife_median), size=1)+
  geom_ribbon(data=m.newdat, aes(x=age,y=ife_median, ymin=m.plo, ymax=m.phi), alpha=0.2)+theme_bw()+
  xlab("Age")+ylab("")+
  theme(axis.text=element_text(size=22),axis.title=element_text(size=22,face="bold"))
ggsave('DIR',width=10, height=8)

#-
# Components separately
#-
m.newdat <- expand.grid(
  age=unique(rr$age), 
  ife_median= 0
) 

m.mm <- model.matrix(terms(ife_age_cont),m.newdat) 
m.newdat$ife_median <- m.mm %*% fixef(ife_age_cont) 
m.pvar1 <- diag(m.mm %*% tcrossprod(vcov(ife_age_cont),m.mm)) 
m.newdat <- data.frame(
  m.newdat
  , m.plo = m.newdat$ife_median-2*sqrt(m.pvar1)
  , m.phi = m.newdat$ife_median+2*sqrt(m.pvar1)

)
m.newdat$ife_median1 <- m.mm[,c(1,2)] %*% fixef(ife_age_cont)[c(1,2)] 
m.newdat$ife_median2 <- m.mm[,c(1,3)] %*% fixef(ife_age_cont)[c(1,3)] 
m.newdat$ife_median3 <- m.mm[,c(1,4)] %*% fixef(ife_age_cont)[c(1,4)] 

ggplot(table2, aes(age,ife_median))+ geom_point(alpha=0.3, size=0.8) + 
  geom_line(data=m.newdat,aes(x=age, y=ife_median), size=1)+
  geom_ribbon(data=m.newdat, aes(x=age,y=ife_median, ymin=m.plo, ymax=m.phi), alpha=0.2)+theme_bw()+
  xlab("Age")+ylab("")+
  theme(axis.text=element_text(size=22),axis.title=element_text(size=22,face="bold"),legend.position="none")+
  geom_line(data=m.newdat,aes(x=age, y=ife_median1, col='1linear'), size=0.8, linetype=2, alpha=0)+
  geom_line(data=m.newdat,aes(x=age, y=ife_median2, col='1quadratic'), size=0.8, linetype=2, alpha=1)+
  geom_line(data=m.newdat,aes(x=age, y=ife_median3, col='cubic'), size=0.8, linetype=2, alpha=0)
ggsave('DIR.pdf',width=10, height=8)


#-
# GENDER
#-
# Modeling 
ife_gender=lmer(ife_median ~ age4*Gender + (1|School), table2)
Anova(ife_gender, type="III") 
emmeans(ife_gender, pairwise ~ Gender, adjust="Bonferroni", type ="response") 
gen.ie=data.frame(emmeans(ife_gender, pairwise ~ Gender,adjust="Bonferroni", type="response")$emmeans)

# Plotting
ggplot(table2, aes(x=Gender, y=ife_median, fill=Gender)) +
  geom_violin(alpha=0.8)+theme_bw()+#geom_boxplot(width=0.08, fill='white')+
  geom_point(position=position_jitter(width=0.1, height=0), size=1, alpha=0.2)+ylab("Inverse efficiency")+xlab("")+
  theme(axis.text=element_text(size=28),axis.title=element_text(size=28,face="bold"))+
  guides(fill=FALSE) +
  geom_point(data=gen.ie, aes(Gender, emmean), shape=15, size=3)+
  geom_errorbar(data=gen.ie, aes(Gender, emmean, ymin = lower.CL, ymax = upper.CL),width = 0.05,size  = 0.5)
ggsave('DIR',width=10, height=8)

#-
# EFFECTS OF AGE AND GENDER ON NUMBER OF TRIALS COMPLETED----------------------------------------------------------------------
#-

tc <- summaryBy(Correct ~ sub + School + age4 + age + Gender, rr, FUN=c(length,mean))
contrasts(tc$age4) = contr.helmert(4)
contrasts(tc$Gender) = contr.helmert(2)

#-
# AGE DISCRETE
#-
# Modeling
trial.comp=lmer(Correct.length ~ age4 + (1|School), tc)
Anova(trial.comp, type="III") 
emmeans(trial.comp, pairwise ~ age4, adjust="Bonferroni") 
emmeans(trial.comp, pairwise ~ age4, adjust="none")

# Tables
cons=as.data.frame(emmeans(trial.comp, pairwise ~ age4, adjust="Bonferroni")$contrasts)
cons$p.uncorrected=round(as.data.frame(emmeans(trial.comp, pairwise ~ age4, adjust="none")$contrasts)$p.value,3)
cons$estimate=round(cons$estimate,2)
cons$SE=round(cons$SE,2)
cons$p.value=round(cons$p.value,3)
cons=cons[,-which(colnames(cons)%in% c("t.ratio", "df"))]
write.csv(cons, "DIR", row.names=F)

# Plots
m.newdat <- expand.grid(
  age4=levels(tc$age4), 
  Correct.length = 0
) 
#contrasts(m.newdat$age4) = contr.helmert(4)
m.mm <- model.matrix(terms(trial.comp),m.newdat) 
m.newdat$Correct.length <- m.mm %*% fixef(trial.comp) 
m.pvar1 <- diag(m.mm %*% tcrossprod(vcov(trial.comp),m.mm)) 
m.newdat <- data.frame(
  m.newdat
  , m.plo = m.newdat$Correct.length-2*sqrt(m.pvar1)
  , m.phi = m.newdat$Correct.length+2*sqrt(m.pvar1)
  #, m.tlo = m.newdat$respC-2*sqrt(m.tvar1)
  #, m.thi = m.newdat$respC+2*sqrt(m.tvar1)
)

ggplot(tc, aes(x=age4, y=Correct.length, fill=age4)) +
  geom_violin(alpha=0.8)+theme_bw()+#geom_boxplot(width=0.08, fill='white')+
  geom_point(position=position_jitter(width=0.1, height=0), size=1, alpha=0.2)+ylab("Number of trials completed")+xlab("Age Group")+
  theme(axis.text.x=element_text(size=18),axis.text.y=element_text(size=22),axis.title=element_text(size=22,face="bold"))+
  scale_fill_discrete(name = "Age group")+guides(fill=FALSE) +
  geom_point(data=m.newdat, aes(age4, Correct.length), shape=15, size=3)+
  geom_errorbar(data=m.newdat, aes(age4, Correct.length, ymin = m.plo, ymax = m.phi),width = 0.05,size  = 0.5)
ggsave('DIR',width=10, height=8)

#-
# AGE CONTINUOUS
#-
# Modeling
trial.comp.cont=lmer(Correct.length ~ poly(age,3) + (1|School), tc)
Anova(trial.comp.cont, type="III") 
summary(trial.comp.cont)

m.newdat <- expand.grid(
  age=unique(tc$age), 
  Correct.length= 0
) 

m.mm <- model.matrix(terms(trial.comp.cont),m.newdat) 
m.newdat$Correct.length <- m.mm %*% fixef(trial.comp.cont) 
m.pvar1 <- diag(m.mm %*% tcrossprod(vcov(trial.comp.cont),m.mm)) 
m.newdat <- data.frame(
  m.newdat
  , m.plo = m.newdat$Correct.length-2*sqrt(m.pvar1)
  , m.phi = m.newdat$Correct.length+2*sqrt(m.pvar1)
)

ggplot(tc, aes(age,Correct.length))+ geom_point(alpha=0.3, size=0.8) + 
  geom_line(data=m.newdat,aes(x=age, y=Correct.length), size=1)+
  geom_ribbon(data=m.newdat, aes(x=age,y=Correct.length, ymin=m.plo, ymax=m.phi), alpha=0.2)+theme_bw()+
  xlab("Age")+ylab("")+
  theme(axis.text=element_text(size=22),axis.title=element_text(size=22,face="bold"))
ggsave('DIR',width=10, height=8)

#-
# Components separately
#-
m.newdat <- expand.grid(
  age=unique(tc$age), 
  Correct.length= 0
) 

m.mm <- model.matrix(terms(trial.comp.cont),m.newdat) 
m.newdat$Correct.length <- m.mm %*% fixef(trial.comp.cont) 
m.pvar1 <- diag(m.mm %*% tcrossprod(vcov(trial.comp.cont),m.mm)) 
m.newdat <- data.frame(
  m.newdat
  , m.plo = m.newdat$Correct.length-2*sqrt(m.pvar1)
  , m.phi = m.newdat$Correct.length+2*sqrt(m.pvar1)
)
m.newdat$length.1 <- m.mm[,c(1,2)] %*% fixef(trial.comp.cont)[c(1,2)] 
m.newdat$length.2 <- m.mm[,c(1,3)] %*% fixef(trial.comp.cont)[c(1,3)] 
m.newdat$length.3 <- m.mm[,c(1,4)] %*% fixef(trial.comp.cont)[c(1,4)] 


ggplot(tc, aes(age,Correct.length))+ geom_point(alpha=0.3, size=0.8) + 
  geom_line(data=m.newdat,aes(x=age, y=Correct.length), size=1)+
  geom_ribbon(data=m.newdat, aes(x=age,y=Correct.length, ymin=m.plo, ymax=m.phi), alpha=0.2)+theme_bw()+
  geom_line(data=m.newdat,aes(x=age, y=length.1, col='1linear'), size=0.8, linetype=2, alpha=0)+
  geom_line(data=m.newdat,aes(x=age, y=length.2, col='2quadratic'), size=0.8, linetype=2)+
  geom_line(data=m.newdat,aes(x=age, y=length.3, col='3cubic'), size=0.8, linetype=2, alpha=0)+theme_bw()+
  xlab("Age")+ylab("")+
  theme(axis.text=element_text(size=22),axis.title=element_text(size=22,face="bold"), legend.position="none")
ggsave('DIR',width=10, height=8)

#-
# GENDER
#-
# Modeling 
tc_gender=lmer(Correct.length ~ age4*Gender + (1|School), tc)
Anova(tc_gender, type="III") 
emmeans(tc_gender, pairwise ~ Gender, adjust="Bonferroni", type ="response") 
gen.ie=data.frame(emmeans(tc_gender, pairwise ~ Gender,adjust="Bonferroni", type="response")$emmeans)

# Plotting
ggplot(tc, aes(x=Gender, y=Correct.length, fill=Gender)) +
  geom_violin(alpha=0.8)+theme_bw()+#geom_boxplot(width=0.08, fill='white')+
  geom_point(position=position_jitter(width=0.1, height=0), size=1, alpha=0.2)+ylab("Number of trials completed")+xlab("")+
  theme(axis.text=element_text(size=28),axis.title=element_text(size=28,face="bold"))+
  guides(fill=FALSE) +
  geom_point(data=gen.ie, aes(Gender, emmean), shape=15, size=3)+
  geom_errorbar(data=gen.ie, aes(Gender, emmean, ymin = lower.CL, ymax = upper.CL),width = 0.05,size  = 0.5)
ggsave('DIR',width=10, height=8)


# PSYCHOMETRICS on reduced dataset (KR-20, p-values, biserial correlations) ----
# Or reduced data set: 
rr.r=read.csv("DIR", sep=",")
missing=rr.r[is.na(rr.r$puzzleSet),c('sub','Trial.Number', 'barcode', 'Relational.Reasoning..Distractor.Strategy','puzzleSet','Relational.Reasoning..Shape.Set')]
length(na.omit(unique(rr.r$barcode))) # How many barcodes is this focusing on? 
sort(unique(rr.r$Trial.Number)) # Which trial.numbers
length(which(rr.r$Correct[rr.r$Reaction.End<250])) # How many decisions with RTs < 250 were correct? 2
# These should be changed to FALSE
rr.r$Correct[rr.r$Correct==T & rr.r$Reaction.End < 250]=FALSE

# How many ppts? 
length(unique(droplevels(na.omit(rr.r$sub))))
# KR-20 (Kuder Richardson formula)
m=dcast(rr.r, sub + age4~Trial.Number, value.var=c("CorrectF"))
round(CronbachAlpha(m[,3:ncol(m)], cond=FALSE, na.rm = TRUE, conf.level=0.95),2)
xtabs(~Trial.Number, rr.r)

for (i in c(1:nrow(trialMap))){
  # The command below merges 2 dfs: the first is a nsub * 2 df, with one column
  # indicating the ppt, and the second indicating the score on all items except
  # the one under examination. The second is a nsub*2 df with one col indicating
  # the ppt and the second whether the item under examination was answered to
  # correctly (true) or not (false)
  t=merge(summaryBy(Correct ~ sub, droplevels(rr.r[rr.r$barcode!=trialMap$barcode[i],]), FUN=mean), 
          rr.r[rr.r$barcode==trialMap$barcode[i],c("Correct", "sub")],all.x=T)
  # Adding a column where correct and incorrect as labaled 1 and 0, respectively
  t$Correct.n=as.numeric(t$Correct) 
  # Computing on version of biserial correlation
  trialMap$biserial[i]=round(biserial.cor(t$Correct.mean, t$Correct, use="complete.obs"),2)
  # Computing simple correlation (which should mathematically be equivalent to the biserial one)
  trialMap$cor[i]=round(cor(t$Correct.mean, t$Correct.n, use="pairwise.complete.obs"),2)
  # Also cross referencing this with the correlation from cor.test, also to get p.value
  trialMap$cor.test[i]=round(cor.test(t$Correct.mean, t$Correct.n, na.action="na.rm")$estimate,2)
  trialMap$cor.test.p[i]=round(cor.test(t$Correct.mean, t$Correct.n, na.action="na.rm")$p.value,3)
  trialMap$n[i]=nrow(t)
}
trialMap=trialMap[order(trialMap$Trial.Number),]
nrow(t)
# P-value info: 
round(mean(trialMap$Correct.mean),2)
round(std.error(trialMap$Correct.mean),2)

# biserial info: 
round(mean(trialMap$cor),2)
round(std.error(trialMap$cor),2)

# Differential item functioning (DIF)
# Age categorical
names=c(levels(m$age4))[-4] 
trialTab.2<- summaryBy(Correct ~ sub, rr, FUN=c(mean))
m.red.2=merge(m, trialTab.2, all.x=T)
difGen.both.2=difGenLogistic(m.red.2[,-c(1,ncol(m.red.2))], group = 'age4', 
                             match=m.red.2$Correct.mean, focal.names = names,
                             p.adjust.method = "BH", type="both")
items.difGen.both.2=difGen.both.2$DIFitems
items.difGen.both.2

# No DIF item detected

# Age continuous
trialTab.2<- summaryBy(Correct ~ sub + age, rr, FUN=c(mean))
m.red.2=merge(m, trialTab.2, all.x=T)
difGen.both.3=difLogistic(m.red.2[,-c(1,2,ncol(m.red.2))], group = 'age', member.type = "cont", 
                          match=m.red.2$Correct.mean, p.adjust.method = "BH", type="both")
items.difGen.both.3=difGen.both.3$DIFitems
items.difGen.both.3

# Gender
trialTab.2<- summaryBy(Correct ~ sub + Gender, rr, FUN=c(mean))
m.red.2=merge(m, trialTab.2, all.x=T)
difGen.both.Gen=difGenLogistic(m.red.2[,-c(1,ncol(m.red.2))], group = 'Gender', 
                               match=m.red.2$Correct.mean, focal.name = 'Male',
                               p.adjust.method = "BH", type="both")
difGen.both.Gen$DIFitems


