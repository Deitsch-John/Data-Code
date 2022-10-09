# Taxa specific results

#Required R packages
library(lme4)
library(tidyverse)
library(gt)
library(MuMIn)
library(nlme)

#Load Datasets
Weather <- read.csv("Weather.csv", header = TRUE) #weather data
Malaise = read.csv("Malaise.csv", header = TRUE) #malaise data
Pans <- read.csv("Pan.Trap.Final.csv", header = TRUE)

# Pan Trap Modeling by Taxa -----------------------------------------------

Pan.Final.Weather <- Pans %>%
  dplyr::inner_join(Weather, by = "study.day")

PanByTaxa <- Pan.Final.Weather %>%
  dplyr::mutate(Parasitoids = (Hymenoptera - Pompilidae - Formicidae + Tachinidae),
         Ichneumonoidea = (Braconidae+Ichneumonidae)) %>%
  dplyr::select(study.day, Group, Plot, Treatment,
         Hours, mean.4day.temp, mean.4day.rain, Total, 
         Hymenoptera, Pompilidae, Formicidae, Parasitoids,
         Ichneumonoidea, Ichneumonidae, Braconidae, 
         Diaprioidea, Chalcidoidea, Coleoptera, Staphylinidae) %>% #select columns of response and predictors
  dplyr::mutate(Total = Total*24/Hours, #scale all taxa groups by 24 hours
         Hymenoptera = Hymenoptera*24/Hours,
         Pompilidae = Pompilidae*24/Hours,
         Formicidae = Formicidae*24/Hours,
         Parasitoids = Parasitoids*24/Hours,
         Ichneumonoidea = Ichneumonoidea*24/Hours,
         Ichneumonidae = Ichneumonidae*24/Hours,
         Braconidae = Braconidae*24/Hours,
         Diaprioidea = Diaprioidea*24/Hours,
         Chalcidoidea = Chalcidoidea*24/Hours,
         Coleoptera = Coleoptera*24/Hours,
         Staphylinidae = Staphylinidae*24/Hours)%>%
  tidyr::pivot_longer(cols = Total:Staphylinidae, "Taxa", values_to = "Abundance")

#function to run models for each taxa, returns a dataframe with model parameters and conf intervals
RunModel_Pan <- function(Taxa_Model){
  taxadata <- dplyr::filter(PanByTaxa, Taxa == Taxa_Model)
  
  fit <- MuMIn::dredge(nlme::gls(scale(log10(Abundance+0.01))~Treatment*scale(study.day)+Group+
                      scale(mean.4day.rain), data = taxadata, na.action = "na.fail", 
                    correlation = nlme::corAR1(, form=~study.day|Plot)))
  
  top_perform <- MuMIn::model.sel(MuMIn::get.models(fit, subset = delta < 2))%>%
    as.data.frame()%>%
    nrow()
  
  coef <- if(top_perform>1){
    as.data.frame(summary(model.avg(fit, subset = delta < 2))$coefmat.subset)
  } else {
    as.data.frame(summary(model.avg(fit))$coefmat.subset)
  }
  
  # CIs <- if(top_perform>1){
  #   as.data.frame(confint(summary(model.avg(fit, subset = delta < 2))))
  # } else {
  #   as.data.frame(confint(summary(model.avg(fit))))
  # }
  
  coef <- as.data.frame(summary(MuMIn::model.avg(fit))$coefmat.subset)
  CIs <- as.data.frame(confint(summary(MuMIn::model.avg(fit))))

  coef2 <- coef %>%
    tibble::rownames_to_column(.,"Parameter")
  
  CIs2 <- CIs %>%
    tibble::rownames_to_column(., "Parameter")%>%
    dplyr::rename(lower = "2.5 %",
           upper = "97.5 %")
  
  info <- dplyr::inner_join(coef2, CIs2, by = "Parameter")%>%
    dplyr::filter(Parameter=="TreatmentIlluminated")%>%
    dplyr::mutate(Taxa=Taxa_Model)
  
  return(info)
  }

TaxaList <- c("Total", "Pompilidae", "Formicidae", "Parasitoids", 
              "Ichneumonoidea", "Diaprioidea", "Chalcidoidea",
              "Coleoptera")

PanParams <- list()
for (i in TaxaList){
  PanParams[[(length(PanParams)+1)]] <- RunModel_Pan(i)
}

ParamsDf <- dplyr::tibble(PanParams)%>%
  tidyr::unnest(PanParams)

#figure 4c
ggplot(data = ParamsDf, aes(x = reorder(Taxa, Estimate), Estimate))+
  geom_hline(yintercept = 0, color = "Red")+
  geom_point(size=5)+
  geom_segment(aes(y=lower, yend=upper,xend=Taxa))+
  coord_flip()+
  labs(x = "Taxa Group", y ="Estimate (95% CI)")+
  labs(x = "Taxa Group", y ="Estimate (95% CI)")+
  theme_classic()+
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11))

#table s4
ParamsDf %>%
  dplyr::filter(Parameter == "TreatmentIlluminated")%>%
  dplyr::select(Taxa, Estimate, `Adjusted SE`:'Pr(>|z|)', lower, upper)%>%
  dplyr::arrange(-Estimate)%>%
  dplyr::mutate(lower = round(lower, digits = 2),
         upper = round(upper, digits = 2),
         Estimate = round(Estimate, digits = 2),
         `Adjusted SE` = round(`Adjusted SE`, digits = 2),
         `z value` = round(`z value`, digits = 2))%>%
  tidyr::unite("CI", lower:upper, remove = TRUE, sep = ", ")%>%
  tidyr::unite("Estimate + SE", Estimate:`Adjusted SE`, remove = TRUE, sep = " + ")%>%
  dplyr::select(Taxa, `Estimate + SE`, CI, `z value`, `Pr(>|z|)`)%>%
  gt::gt()%>%
  gt::fmt_number(
    columns = "Pr(>|z|)",
    decimals = 3
  )

# Malaise Trap Modeling by Taxa -----------------------------------------------

MalaiseByTaxa <- Malaise %>%
  dplyr::mutate(Predators = Spiders + Pompilidae + Formicidae + Asilidae,
         PredPtoid = Parasitoids + Predators,
         PredProp = Predators/Total,
         PredPtoidProp = PredPtoid/Total,
         Ichneumonoidea = Ichneumonidae + Braconidae,
         MothProp = Lepidoptera/Total) %>%
  dplyr::mutate(Date = mdy(Date),
         yday = yday(Date),
         study.day = yday - 151,
         day.postlights = study.day - 13) %>%
  dplyr::inner_join(Weather, by = "study.day")%>%
  dplyr::select(day.postlights, Group, Plot, Treatment,
         Hours, mean.4day.temp, mean.4day.rain, Total,
         Predators, Parasitoids, PredPtoid, PredProp, PredPtoidProp, MothProp,
         Lepidoptera, Plecoptera, Coleoptera, Spiders, Diptera, Tipuloidea,
         Hymenoptera, Ichneumonoidea, Ichneumonidae, Braconidae,
         Diaprioidea, Formicidae, Pompilidae) %>%
  dplyr::mutate(Total = Total*24/Hours,
         Predators = Predators*24/Hours,
         Parasitoids = Parasitoids*24/Hours,
         PredPtoid = PredPtoid*24/Hours,
         Lepidoptera = Lepidoptera*24/Hours,
         Plecoptera = Plecoptera*24/Hours,
         Coleoptera = Coleoptera*24/Hours,
         Spiders = Spiders*24/Hours,
         Diptera = Diptera*24/Hours,
         Tipuloidea = Tipuloidea*24/Hours,
         Hymenoptera = Hymenoptera*24/Hours,
         Ichneumonoidea = Ichneumonoidea*24/Hours,
         Ichneumonidae = Ichneumonidae*24/Hours,
         Braconidae = Braconidae*24/Hours,
         Diaprioidea = Diaprioidea*24/Hours,
         Pompilidae = Pompilidae*24/Hours,
         Formicidae = Formicidae*24/Hours) %>%
  dplyr::filter(day.postlights >= 0)%>%
  tidyrr::pivot_longer(Total:Pompilidae, names_to = "Taxa", values_to = "Abundance")

#function to run models for each taxa, returns a dataframe with model parameters and conf intervals
RunModel_Malaise <- function(Taxa_Model){
  taxadata <- dplyr::filter(MalaiseByTaxa, Taxa == Taxa_Model)
  
  fit <- MuMIn::dredge(nlme::gls(scale(log10(Abundance+0.01))~Treatment*scale(day.postlights)+Group+
                      scale(mean.4day.rain)+scale(mean.4day.temp), 
                    data = taxadata, na.action = "na.fail", correlation = 
                      nlme::corAR1(, form=~day.postlights|Plot)))
  
  top_perform <- MuMIn::model.sel(MuMIn::get.models(fit, subset = delta < 2))%>%
    as.data.frame()%>%
    nrow()
  
  # coef <- if(top_perform>1){
  #   as.data.frame(summary(model.avg(fit, subset = delta < 2))$coefmat.subset)
  # } else {
  #   as.data.frame(summary(model.avg(fit))$coefmat.subset)
  # }
  # 
  # CIs <- if(top_perform>1){
  #   as.data.frame(confint(summary(model.avg(fit, subset = delta < 2))))
  # } else {
  #   as.data.frame(confint(summary(model.avg(fit))))
  # }
  
  coef <- as.data.frame(summary(MuMIn::model.avg(fit))$coefmat.subset)
  CIs <- as.data.frame(confint(summary(MuMIn::model.avg(fit))))
  
  coef2 <- coef %>%
    tibble::rownames_to_column(.,"Parameter")
  
  CIs2 <- CIs %>%
    tibble::rownames_to_column(., "Parameter")%>%
    dplyr::rename(lower = "2.5 %",
           upper = "97.5 %")
  
  info <- dplyr::inner_join(coef2, CIs2, by = "Parameter")%>%
    dplyr::filter(Parameter=="TreatmentIlluminated")%>%
    dplyr::mutate(Taxa=Taxa_Model)
  
  return(info)
}

TaxaList_Malaise <- c("Total", "Pompilidae", "Formicidae", 
              "Ichneumonoidea", "Diaprioidea", "Coleoptera",
              "Diptera", "Tipuloidea", "Spiders", "Lepidoptera")

MalaiseParams <- list()
for (i in TaxaList_Malaise){
  MalaiseParams[[(length(MalaiseParams)+1)]] <- RunModel_Malaise(i)
}

ParamsDf_Mal <- dplyr::tibble(MalaiseParams)%>%
  tidyr::unnest(MalaiseParams)

#figure 4f
ggplot(data = ParamsDf_Mal, aes(x = reorder(Taxa, Estimate), Estimate))+
  geom_hline(yintercept = 0, color = "Red")+
  geom_point(size=5)+
  geom_segment(aes(y=lower, yend=upper,xend=Taxa))+
  coord_flip()+
  labs(x = "Taxa Group", y ="Estimate (95% CI)")+
  theme_classic()+
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11))

#table s7
ParamsDf_Mal %>%
  dplyr::filter(Parameter == "TreatmentIlluminated")%>%
  dplyr::select(Taxa, Estimate, `Adjusted SE`:'Pr(>|z|)', lower, upper)%>%
  dplyr::arrange(-Estimate)%>%
  dplyr::mutate(lower = round(lower, digits = 2),
                upper = round(upper, digits = 2),
                Estimate = round(Estimate, digits = 2),
                `Adjusted SE` = round(`Adjusted SE`, digits = 2),
                `z value` = round(`z value`, digits = 2))%>%
  tidyr::unite("CI", lower:upper, remove = TRUE, sep = ", ")%>%
  tidyr::unite("Estimate + SE", Estimate:`Adjusted SE`, remove = TRUE, sep = " + ")%>%
  dplyr::select(Taxa, `Estimate + SE`, CI, `z value`, `Pr(>|z|)`)%>%
  gt::gt()%>%
  gt::fmt_number(
    columns = "Pr(>|z|)",
    decimals = 3
  )

