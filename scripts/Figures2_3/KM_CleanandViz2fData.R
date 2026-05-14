setwd("~/CM_515/FinalProj/data")

library(tidyverse)
library(readxl)
library(ggpubr)
#install.packages("egg")
library(egg)
library(ggthemes)

transact <- read_xlsx("D_SourceDataFigs.xlsx", sheet = "Fig. 2f")

#Clean data
transact_clean <- transact%>%
  mutate(effector = gsub("GFP", "Ctrl",effector))%>%#Rename column contents to how they should appear on plot
  mutate(effector = gsub("SlSSP2", "SSP2F169", effector))%>%
  mutate(effector = gsub("SpSSP2", "SSP2S169", effector))%>%
  mutate(promoter = gsub("FUL1", "SlFUL1", promoter))%>%
  mutate(promoter = gsub("FUL2", "SlFUL2", promoter))%>%
  mutate(`fLUC/rLUC` = as.numeric(`fLUC/rLUC`))%>%#Convert column to numeric and then filter out NAs
  filter(!is.na(`fLUC/rLUC`))%>%
  mutate(effector = as.factor(effector))%>%#Convert effector column to factor data and relevel factors to order they need to appear on plot
  mutate(effector = fct_relevel(effector, "SSP2F169","SSP2S169", "SSP", "Ctrl"))

#Plot
#Add data points on top of boxplots with layered geometries: geom_boxplot and geom_jitter
#Use coordflip() to flip the plot on its side (so y axis is along the bottom and x axis is along the side)

#CAN PLOT TWO WAYS

#PLOT VIA FACET WRAP
figure2f_facet <- transact_clean%>%
  ggplot(aes(x = effector, y = `fLUC/rLUC`))+
  facet_wrap(vars(promoter), labeller = as_labeller(c("pMC" = "pMC::fLUC", "SlFUL1" = "SlFUL1::fLUC", "SlFUL2" = "SlFUL2::fLUC")), ncol = 1)+#Set labels of each plot
  #and direct the function to arrange the plots in 1 column
  geom_boxplot(coef = NULL)+
  geom_jitter(color = "steelblue1", size = 1)+
  scale_y_continuous(limits = c(0, 2))+
  theme_few()+
  coord_flip()+
  labs(x = "", y = "Luciferase ratio (fLUC/rLUC)")+
  theme(strip.text = element_text(hjust = 0))

#PLOT VIA GGARRANGE
#Create  separate dfs for each promoter
#Create boxplots for each promoter
  #Use theme(plot.margin) to make sure plots space out well enough when grouped with ggarrange()
  #Remove x axis labels from all 3 plots, remove y axis labels from all except the bottom

transact_pMC <- transact_clean%>%
  filter(promoter == "pMC")

pMCplot<- transact_pMC%>%
  ggplot(aes(x = effector, y = `fLUC/rLUC`))+
  geom_boxplot(coef = NULL)+
  geom_jitter(color = "steelblue1", size = 1)+
  scale_y_continuous(limits = c(0, 2))+
  theme_classic()+
  coord_flip()+
  theme(plot.margin = margin(t = 20, r = 20, b = 20, l = 20))+
  labs(x = "", y = "")
  
transact_SlFUL1 <- transact_clean%>%
  filter(promoter == "SlFUL1")

SlFUL1plot <- transact_SlFUL1%>%
  ggplot(aes(x = effector, y = `fLUC/rLUC`))+
  geom_boxplot(coef = NULL)+
  geom_jitter(color = "steelblue1", size = 1)+
  scale_y_continuous(limits = c(0, 2))+
  theme_classic()+
  coord_flip()+
  theme(plot.margin = margin(t = 20, r = 20, b = 20, l = 20))+
  labs(x = "", y = "")

transact_SlFUL2 <- transact_clean%>%
  filter(promoter == "SlFUL2")

SlFUL2plot<- transact_SlFUL2%>%
  ggplot(aes(x = effector, y = `fLUC/rLUC`))+
  geom_boxplot(coef = NULL)+
  geom_jitter(color = "steelblue1", size = 1)+
  scale_y_continuous(limits = c(0, 2))+
  theme_classic()+
  coord_flip()+
  theme(plot.margin = margin(t = 20, r = 20, b = 20, l = 20))+
  labs(y = "Luciferase ratio (fLUC/rLUC)", x = "")


#Arrange 3 plots into one figure
figure2f_arrange <- ggpubr::ggarrange(pMCplot, SlFUL1plot, SlFUL2plot,
                      labels = c("pMC::fLUC", "pSlFUL1::fLUC", "pSlFUL2::fLUC"),
                      ncol = 1, nrow = 3,
                      font.label = list(size = 10))

figure2f_facet
figure2f_arrange