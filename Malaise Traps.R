# Malaise Trap Data Analysis

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
Malaise <- read.csv("Malaise_Raw.csv", header = TRUE)
Weather = read.csv("Weather.csv", header = TRUE)

MalaiseByTaxa <- Malaise %>%
  mutate(Predators = Spiders + Pompilidae + Formicidae + Asilidae,
         PredPtoid = Parasitoids + Predators,
         PredProp = Predators/Total,
         PredPtoidProp = PredPtoid/Total,
         Ichneumonoidea = Ichneumonidae + Braconidae,
         MothProp = Lepidoptera/Total) %>%
  mutate(Date = mdy(Date),
         yday = yday(Date),
         study.day = yday - 151,
         day.postlights = study.day - 13) %>%
  inner_join(Weather, by = "study.day")%>%
  mutate(Total = Total*24/Hours,
         Predators = Predators*24/Hours,
         Parasitoids = Parasitoids*24/Hours)%>%
  select(Sample, Group, Pair, Plot, Treatment, Road, 
         Date, Hours, Predators, Parasitoids, yday, study.day,
         mean.4day.temp, mean.4day.rain, day.postlights)%>%
  filter(day.postlights > 0)

# Initial Data expl. ---------------------------------------------------------------
hist(log(MalaiseByTaxa$Predators))
hist(log(MalaiseByTaxa$Parasitoids))

plot(density(log(MalaiseByTaxa$Predators)))
plot(density(log(MalaiseByTaxa$Parasitoids)))

summary(aov(log10(Predators + 0.1)~Group, data = MalaiseByTaxa))
summary(aov(log10(Parasitoids + 0.1)~Group, data = MalaiseByTaxa))
summary(aov(log10(Predators + 0.1)~Pair, data = MalaiseByTaxa))
summary(aov(log10(Parasitoids + 0.1)~Pair, data = MalaiseByTaxa)) #significant
summary(aov(log10(Predators + 0.1)~Road, data = MalaiseByTaxa))
summary(aov(log10(Parasitoids + 0.1)~Road, data = MalaiseByTaxa)) #significant
summary(aov(log10(Predators + 0.1)~Sample, data = MalaiseByTaxa)) #significant
summary(aov(log10(Parasitoids + 0.1)~Sample, data = MalaiseByTaxa))#significant

# linear mod --------------------------------------------------------------
malaise.lm <- lm(log10(Predators+0.1) ~ Treatment*scale(day.postlights) +
                   scale(mean.4day.temp) + scale(mean.4day.rain), 
                 data = MalaiseByTaxa, na.action = "na.fail")
#check for correlation among predictors
vif(malaise.lm) #vif < 3

#checking for autocorrelation
durbinWatsonTest(malaise.lm) #no autocorrelation
acf(malaise.lm$residuals) #no autocorrelation!
resid_panel(malaise.lm)


# gls malaise pred --------------------------------------------------------

malaise.gls.pred <- gls(log10(Predators+0.1) ~ Treatment*scale(day.postlights) + 
                          scale(mean.4day.temp) + scale(mean.4day.rain), 
                        data=MalaiseByTaxa, 
                        correlation = corAR1(, form=~day.postlights|Plot))

vif(malaise.gls.pred)

malaise.pred <- dredge(malaise.gls.pred)
ModelTable.pred <- model.sel(get.models(malaise.pred, subset = TRUE))

#table s5
ModelTable.pred %>%
  as.data.frame()%>%
  mutate(modelno = rownames(ModelTable.pred))%>%
  select(-Treatment,-`scale(mean.4day.rain)`,
         -`scale(mean.4day.temp)`,-`(Intercept)`,
         -`scale(day.postlights)`, -`scale(day.postlights):Treatment`)%>%
  # filter(delta < 10)%>%
  gt(rowname_col = "modelno")%>%
  fmt_number(
    columns = everything(),
    decimals = 2)%>%
  cols_width(everything()~px(100))

modavgd.pred <- model.avg(ModelTable.pred, subset = delta < 2) #1 best model
modavgd.summary.pred <- summary(modavgd.pred)
modavgd.coeff.pred <- as.data.frame(modavgd.summary.pred$coefmat.subset)
modavgd.CI.pred <- as.data.frame(confint(modavgd.summary.pred)) #confidence int
sw(malaise.pred)

modavdg.details.pred <- bind_cols(modavgd.coeff.pred, modavgd.CI.pred)

#Table for Manuscript; view in html and copy to word doc
modavdg.details.pred %>%
  mutate(name = rownames(modavdg.details.pred))%>%
  rename(lower = "2.5 %",
         upper = "97.5 %")%>%
  mutate(lower = round(lower, digits = 4),
         upper = round(upper, digits = 4))%>%
  unite("CI", lower:upper, remove = TRUE, sep = ", ")%>%
  gt(rowname_col = "name") %>%
  fmt_number(
    columns = c("Estimate", "Std. Error", "Adjusted SE",
                "z value", "Pr(>|z|)"),
    decimals = 4)%>%
  tab_stubhead(label = "Parameter")%>%
  cols_width(c("Estimate", "Std. Error", "Adjusted SE",
               "z value", "Pr(>|z|)")~px(100))


hist(malaise.gls.pred$residuals)

# gls malaise parasitoids -------------------------------------------------

malaise.gls.par <- gls(log10(Parasitoids+0.1) ~ Treatment*scale(day.postlights) + 
                         scale(mean.4day.temp) + scale(mean.4day.rain) + Group + Road, 
                       data=MalaiseByTaxa, 
                       correlation = corAR1(, form=~day.postlights|Plot))

vif(malaise.gls.par)

malaise.par <- dredge(malaise.gls.par)
ModelTable.par <- model.sel(get.models(malaise.par, subset = TRUE))

#table s6
ModelTable.par %>%
  as.data.frame()%>%
  mutate(modelno = rownames(ModelTable.par))%>%
  select(-Treatment,-`scale(mean.4day.rain)`,
         -`scale(mean.4day.temp)`,-`(Intercept)`,
         -`scale(day.postlights)`, -`scale(day.postlights):Treatment`,
         -Road, -Group)%>%
  # filter(delta < 2)%>%
  gt(rowname_col = "modelno")%>%
  fmt_number(
    columns = everything(),
    decimals = 2)%>%
  cols_width(everything()~px(100))

#finding average parameters
modavgd.par <- model.avg(ModelTable.par, subset = delta < 2) #1 best model
modavgd.summary.par <- summary(modavgd.par)
modavgd.coeff.par <- as.data.frame(modavgd.summary.par$coefmat.subset)
modavgd.CI.par <- as.data.frame(confint(modavgd.summary.par)) #confidence int
sw(malaise.par)

modavdg.details.par <- bind_cols(modavgd.coeff.par, modavgd.CI.par)

#table 3
modavdg.details.par %>%
  mutate(name = rownames(modavdg.details.par))%>%
  rename(lower = "2.5 %",
         upper = "97.5 %")%>%
  mutate(lower = round(lower, digits = 4),
         upper = round(upper, digits = 4))%>%
  unite("CI", lower:upper, remove = TRUE, sep = ", ")%>%
  gt(rowname_col = "name") %>%
  fmt_number(
    columns = c("Estimate", "Std. Error", "Adjusted SE",
                "z value", "Pr(>|z|)"),
    decimals = 4)%>%
  tab_stubhead(label = "Parameter")%>%
  cols_width(c("Estimate", "Std. Error", "Adjusted SE",
               "z value", "Pr(>|z|)")~px(100))


#residuals normality
hist(malaise.gls.par$residuals)

# figures -----------------------------------------------------------------

#figure 4D
ggplot(MalaiseByTaxa, aes(Treatment, log10(Predators + 0.1)))+
  geom_boxplot(aes(fill = Treatment))+
  geom_signif(comparisons = list(c("Illuminated", "Control")),
              map_signif_level = TRUE)+
  geom_jitter(width = 0.06, height = 0, alpha = 0.7)+
  scale_y_continuous(limits = c(-1.5,1.5), breaks = seq(-1.5,1.5, by = 1))+
  theme_minimal()+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = "none")+
  labs(y = "log (Predators / 24 hours)",
       x = "Treatment")

#figure 4E
ggplot(MalaiseByTaxa, aes(Treatment, log10(Parasitoids + 0.1)))+
  geom_boxplot(aes(fill = Treatment))+
  geom_signif(comparisons = list(c("Illuminated", "Control")),
              map_signif_level = TRUE)+
  geom_jitter(width = 0.06, height = 0, alpha = 0.7)+
  scale_y_continuous(limits = c(-1.5,1.5), breaks = seq(-1.5,1.5, by = 1))+
  theme_minimal()+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = "none")+
  labs(y = "log (Parasitoids / 24 hours)",
       x = "Treatment")
