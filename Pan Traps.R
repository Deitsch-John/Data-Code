# Pan Trap Data Analysis

# Required R Packages
library(tidyverse)
library(lme4)
library(MuMIn)
library(gt)
library(ggsignif)
library(car)
library(ggResidpanel)
library(lubridate)
library(lawstat)
library(nlme)
library(rstatix)


# Data Prep ---------------------------------------------------------------
Pans <- read.csv("Pan_Raw.csv", header = TRUE) #pan trap data
Weather <- read.csv("Weather.csv", header = TRUE) #weather data

Pan.Final.Weather <- Pans %>%
  mutate(Parasitoids = Ichneumonidae + Braconidae + Chalcidoidea + 
           Diaprioidea + OtherParasitoids + Tachinidae,
         Predators = Pompilidae + Formicidae + Asilidae + Staphylinidae + Carabidae) %>%
  mutate(Parasitoids = (Parasitoids*24/Hours),
         Predators = (Predators*24/Hours))%>%
  inner_join(Weather, by = "study.day")%>%
  mutate(log.parasitoids = log10(Parasitoids + 0.01),
         log.predators = log10(Predators + 0.01))


# Initial data expl. -----------------------------------------------

#Check normality of pan trap abundance 
shapiro.test(Pan.Final.Weather$log.Total) #normal

#Did abundance vary between N and S blocks?
summary(aov(log.parasitoids~Group, data = Pan.Final.Weather))
summary(aov(log.predators~Group, data = Pan.Final.Weather))

#yes

#Did abundance vary between pairs?

summary(aov(log.parasitoids~Group+Pair, data = Pan.Final.Weather))
summary(aov(log.predators~Group+Pair, data = Pan.Final.Weather))

#no

#Did abundance vary significantly between 'road' plots and 'off-road' plots?

summary(aov(log.parasitoids~Road, data = Pan.Final.Weather))
summary(aov(log.predators~Road, data = Pan.Final.Weather))

#no

#repeated measures ANOVA on pred rate with respect to time
Pan.Final.Weather%>%
  anova_test(dv = log.parasitoids, wid = Plot, within = Sample) %>%
  get_anova_table()%>%
  gt()

Pan.Final.Weather%>%
  anova_test(dv = log.predators, wid = Plot, within = Sample) %>%
  get_anova_table()%>%
  gt()

#decreased over time
# linear mod --------------------------------------------------------------

#define global model
global.pan.pr <- lm(log.predators ~ Group + Treatment + scale(study.day) +
                      scale(mean.4day.temp) + scale(mean.4day.rain), 
                    data = Pan.Final.Weather, na.action = "na.fail")
#check for correlation among predictors
vif(global.pan.pr) #vif > 3

#remove study.day
global.pan.pr <- lm(log.parasitoids~Group+Treatment+scale(study.day)+
                      scale(mean.4day.rain), data = Pan.Final.Weather, na.action = "na.fail")

#check for correlation among predictors
vif(global.pan.pr) #good to go

#checking for autocorrelation
durbinWatsonTest(global.pan.pr) #p < 0.05
acf(global.pan.pr$residuals) #lag 1 autocorrelation


# gls pan parasitoid --------------------------------------------------------------

#models with linear distribution have auto-correlated residuals
#gls model with ar1 symmetry fit to try to over come this issue

#auto-regressive symmetry global model
gls.global.pa <- gls(log.parasitoids ~ Group + Treatment*scale(study.day) + scale(mean.4day.temp) + scale(mean.4day.rain), data=Pan.Final.Weather, correlation = corAR1(, form=~study.day|Plot))

vif(gls.global.pa) #a little high

#remove temperature

gls.global.pa <- gls(log.parasitoids ~ Group + Treatment*scale(study.day)+ scale(mean.4day.rain), data=Pan.Final.Weather, correlation = corAR1(, form=~study.day|Plot))

vif(gls.global.pa) #we good

gls.all.pa <- dredge(gls.global.pa)
gls.Table <- model.sel(get.models(gls.all.pa, subset = TRUE))

#Table S3
gls.Table %>%
  as.data.frame()%>%
  mutate(modelno = rownames(gls.Table))%>%
  select(-Group,-Treatment,-`scale(mean.4day.rain)`,
         -`(Intercept)`, -`scale(study.day)`, -`scale(study.day):Treatment`)%>%
  gt(rowname_col = "modelno")%>%
  fmt_number(
    columns = everything(),
    decimals = 2)%>%
  cols_width(everything()~px(120))%>%
  tab_stubhead(label = "Model")%>%
  tab_source_note(source_note = "text")%>%
  tab_header(title = "sample")

#finding average parameters
modavdg.pan.gls.pa <- model.sel(get.models(gls.Table, subset = 1))

gls.top.pa <- gls(log.parasitoids ~ Group + Treatment, data=Pan.Final.Weather, correlation = corAR1(, form=~study.day|Plot))

topmodel.summ <- summary(gls.top.pa)
topcoeff <- as.data.frame(topmodel.summ$tTable)
topci <- as.data.frame(confint(gls.top.pa))

topdetails <- bind_cols(topcoeff, topci)

#table 2
topdetails %>%
  mutate(name = rownames(topdetails))%>%
  rename(lower = "2.5 %",
         upper = "97.5 %",
         Estimate="Value",
         `Std. Error`="Std.Error",
         `z value`="t-value",
         `Pr(>|z|)`="p-value")%>%
  mutate(lower = round(lower, digits = 2),
         upper = round(upper, digits = 2))%>%
  unite("CI", lower:upper, remove = TRUE, sep = ", ")%>%
  gt(rowname_col = "name") %>%
  fmt_number(
    columns = c("Estimate", "Std. Error",
                "z value", "Pr(>|z|)"),
    decimals = 3)%>%
  tab_stubhead(label = "Parameter")%>%
  cols_width(c("Estimate", "Std. Error",
               "z value", "Pr(>|z|)", "CI")~px(100))

#residuals normality
hist(gls.top.pa$residuals)

# gls pan predator --------------------------------------------------------------

gls.global.pr <- gls(log.predators ~ Group + Treatment*scale(study.day) + scale(mean.4day.temp) + scale(mean.4day.rain), data=Pan.Final.Weather, correlation = corAR1(, form=~study.day|Plot))

vif(gls.global.pr) #a little high

#remove temperature

gls.global.pr <- gls(log.predators ~ Group + Treatment*scale(study.day)+ scale(mean.4day.rain), data=Pan.Final.Weather, correlation = corAR1(, form=~study.day|Plot))

vif(gls.global.pr) #we good

#creating candidate set of models
gls.all.pr <- dredge(gls.global.pr)
gls.Table.pr <- model.sel(get.models(gls.all.pr, subset = TRUE))

#table S2
gls.Table.pr %>%
  as.data.frame()%>%
  mutate(modelno = rownames(gls.Table.pr))%>%
  select(-Group,-Treatment,-`scale(mean.4day.rain)`,
         -`(Intercept)`, -`scale(study.day)`, -`scale(study.day):Treatment`)%>%
  filter(delta < 4)%>%
  gt(rowname_col = "modelno")%>%
  fmt_number(
    columns = everything(),
    decimals = 2)%>%
  cols_width(everything()~px(120))%>%
  tab_stubhead(label = "Model")%>%
  tab_source_note(source_note = "text")%>%
  tab_header(title = "sample")

modavgd.pan.gls.pr <- model.avg(gls.Table.pr, subset = delta < 2)
modavgd.summary.pan.gls.pr <- summary(modavgd.pan.gls.pr)
modavgd.coeff.pan.gls.pr <- as.data.frame(modavgd.summary.pan.gls.pr$coefmat.subset)
modavgd.CI.pan.gls.pr <- as.data.frame(confint(modavgd.summary.pan.gls.pr))
sw(gls.all.pr)

modavdg.details.pan.gls.pr <- bind_cols(modavgd.coeff.pan.gls.pr, modavgd.CI.pan.gls.pr)

#table 2
modavdg.details.pan.gls.pr %>%
  mutate(name = rownames(modavdg.details.pan.gls.pr))%>%
  rename(lower = "2.5 %",
         upper = "97.5 %")%>%
  mutate(lower = round(lower, digits = 2),
         upper = round(upper, digits = 2))%>%
  unite("CI", lower:upper, remove = TRUE, sep = ", ")%>%
  gt(rowname_col = "name") %>%
  fmt_number(
    columns = c("Estimate", "Std. Error", "Adjusted SE",
                "z value", "Pr(>|z|)"),
    decimals = 3)%>%
  tab_stubhead(label = "Parameter")%>%
  cols_width(c("Estimate", "Std. Error", "Adjusted SE",
               "z value", "Pr(>|z|)", "CI")~px(100))

#residuals normality
hist(gls.global.pr$residuals)
# figures -----------------------------------------------------------------

#figure 4a
ggplot(Pan.Final.Weather, aes(Treatment, log.predators))+
  geom_boxplot(aes(fill = Treatment))+
  geom_signif(comparisons = list(c("Illuminated", "Control")),
              map_signif_level = TRUE)+
  geom_jitter(width = 0.06, height = 0, alpha = 0.7)+
  scale_y_continuous(limits = c(-2.5,2.5), breaks = seq(-2.5,2.5, by = 1))+
  theme_minimal()+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = "none")+
  labs(y = "log (Predator Abundance)",
       x = "Treatment")

#figure 4b
ggplot(Pan.Final.Weather, aes(Treatment, log.parasitoids))+
  geom_boxplot(aes(fill = Treatment))+
  geom_signif(comparisons = list(c("Illuminated", "Control")),
              map_signif_level = TRUE)+
  geom_jitter(width = 0.06, height = 0, alpha = 0.7)+
  scale_y_continuous(limits = c(-2.5,2.5), breaks = seq(-2.5,2.5, by = 1))+
  theme_minimal()+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = "none")+
  labs(y = "log (Parasitoid Abundance)",
       x = "Treatment")
