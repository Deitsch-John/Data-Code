# Predation Rates

#Required R packages
library(tidyverse)
library(lme4)
library(MuMIn)
library(gt)
library(ggsignif)
library(car)
library(ggResidpanel)
library(nlme)
library(rstatix)

# Data prep -----------------------------------------------

models <- read.csv("Caterpillars_Raw.csv", header = TRUE)

PredRates <- models %>%
  mutate(score = case_when(
    Marked.Arthropod=="Yes"~1,
    Marked.Arthropod=="No"~0
  ))%>%
  group_by(Period, Group, Pair, Plot, Treatment)%>%
  summarize(sample_size = n(),
            pred = sum(score))%>%
  mutate(Prop.marked.arthropod = pred/sample_size)%>%
  mutate(Road = case_when(
    Plot%in%c("Control II", "ALAN III", "ALAN VI", "Control VII")~"No",
    Plot%in%c("Control IV", "ALAN I", "ALAN VIII", "Control V")~"Yes"
  ),
  Moon = case_when(
    Period%in%c("2","4")~"Full",
    Period%in%c("1","3")~"New"
  ),
  Year.Day = case_when(
    Period==1~166,
    Period==2~179,
    Period==3~194,
    Period==4~208,
  ),
  ) %>%
  mutate(Study.Day = Year.Day-151)

ALAN <- PredRates %>%
  ungroup()%>%
  select(Period, Pair, Treatment, Prop.marked.arthropod)%>%
  filter(Treatment=="Illuminated")%>%
  pivot_wider(names_from = Treatment,
              values_from = Prop.marked.arthropod)%>%
  mutate(identifier = str_c(Pair, Period, sep = " "))

Control <- PredRates %>%
  ungroup()%>%
  select(Period, Pair, Treatment, Prop.marked.arthropod)%>%
  filter(Treatment=="Control")%>%
  pivot_wider(names_from = Treatment,
              values_from = Prop.marked.arthropod)%>%
  mutate(identifier = str_c(Pair, Period, sep = " "))%>%
  select(-Pair, -Period)

Paired <- inner_join(ALAN, Control, by = "identifier")%>%
  mutate(Ratio = Illuminated/Control)


# Initial data expl. -----------------------------------------------

shapiro.test(PredRates$Prop.marked.arthropod) 
plot(density(PredRates$Prop.marked.arthropod))

aov.group <- aov(Prop.marked.arthropod~Group, PredRates)
summary(aov.group)
print(TukeyHSD(aov.group, which = "Group"))

aov.pair <- aov(Prop.marked.arthropod~Pair, PredRates)
summary(aov.pair)
print(TukeyHSD(aov.pair, which = "Pair"))

aov.road <- aov(Prop.marked.arthropod~Road, PredRates)
summary(aov.road)
print(TukeyHSD(aov.road, which = "Road"))

PredRates %>%
  anova_test(dv = Prop.marked.arthropod, wid = Plot, within = Period) %>%
  get_anova_table()%>%
  gt()

shapiro.test(Paired$Ratio) #normal
shapiro.test(log10(Paired$Ratio)) #more normal

t.test(Paired$Ratio, mu = 0, alternative = "two.sided")

Paired %>%
  anova_test(dv = Ratio, wid = Pair, within = Period) %>%
  get_anova_table()%>%
  gt()

# Models -----------------------------------------------

PlasticineGlobal = lm(Prop.marked.arthropod~
                        Period + Treatment, data = PredRates, 
                      na.action = "na.fail")

AllModels <- dredge(PlasticineGlobal)
ModelTable <- model.sel(get.models(AllModels, subset = NA))

#Table S1
ModelTable %>%
  as.data.frame()%>%
  mutate(modelno = rownames(ModelTable))%>%
  select(-Period,-Treatment,`(Intercept)`)%>%
  gt(rowname_col = "modelno")%>%
  fmt_number(
    columns = everything(),
    decimals = 2)%>%
  cols_width(everything()~px(120))

topmodel.summ <- summary(PlasticineGlobal)
topcoeff <- as.data.frame(topmodel.summ$coefficients)
topci <- as.data.frame(confint(PlasticineGlobal))

topdetails <- bind_cols(topcoeff, topci)

#Table 1
topdetails %>%
  mutate(name = rownames(topdetails))%>%
  rename(lower = "2.5 %",
         upper = "97.5 %")%>%
  mutate(lower = round(lower, digits = 2),
         upper = round(upper, digits = 2))%>%
  unite("CI", lower:upper, remove = TRUE, sep = ", ")%>%
  gt(rowname_col = "name")%>%
  fmt_number(
    columns = c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
    decimals = 4)

# model eval  -----------------------------------------------

durbinWatsonTest(PlasticineGlobal) #no autocorrelation
acf(PlasticineGlobal$residuals) #autocorrelation plot

#no autocorrelation detected, a few close lines on the ACF plot

#residuals homogeneity and normality
ggResidpanel::resid_panel(PlasticineGlobal) #no alarming pattern  
shapiro.test(PlasticineGlobal$residuals) #normal


# figures -----------------------------------------------

#figure 3a
PredRates%>%
  ggplot(aes(Treatment, Prop.marked.arthropod))+
  geom_boxplot(aes(fill = Treatment))+
  geom_signif(comparisons = list(c("Illuminated", "Control")),
              map_signif_level = TRUE)+
  geom_jitter(width = 0.06, height = 0, alpha = 0.8)+
  scale_y_continuous(limits = c(0.15,0.9), breaks = seq(0.15, 0.9, by = 0.15))+
  theme_minimal()+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))+
  labs(y = "Predation Rate from Arthropods",
       x = "Treatment")

#figure 3b
PredRates %>%
  group_by(Period, Treatment)%>%
  summarise(prop.marked.mean = mean(Prop.marked.arthropod),
            sd.prop = sd(Prop.marked.arthropod),
            sample.size = n())%>%
  mutate(stand.e = sd.prop/sqrt(sample.size))%>%
  rename("Mean Predation Rate" = prop.marked.mean,
         "SD" = sd.prop,
         "Sample Size" = sample.size,
         "Standard Error" = stand.e)%>%
  ungroup()%>%
  ggplot(aes(Period, `Mean Predation Rate`, 
             shape = Treatment, color = Treatment))+
  geom_point(size = 4)+
  geom_line()+
  geom_errorbar(aes(ymin = `Mean Predation Rate`- `Standard Error`,
                    ymax = `Mean Predation Rate`+ `Standard Error`),
                width = 0.2, size = 0.7)+
  theme_minimal()+
  scale_y_continuous(limits = c(0.3, 0.7))+
  labs(y = "Predation Rate (Mean +/- SE)",
       x = "Sampling Period")+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))+
  scale_y_continuous(limits = c(0.15,0.90), breaks = seq(0.15,0.9, by =0.15))

