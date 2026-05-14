##############
##############

setwd("~/CM_515/FinalProj/data")

library(readxl)
library(tidyverse)
library(dplyr)
library(plyr)
library(ggpubr)
library(ggthemes)
#install.packages("ggh4x")
library(ggh4x)
library(S4Vectors)
#This script requires S4Vectors because when I wrote this script I must have had S4Vectors loaded from another script,
#and when I was using the rename() function, the dplyr syntax wasn’t working but the reverse order (S4Vectors syntax) did work.
#As a result, I left the S4Vectors syntax in the script.
#When I ran the script again later during a different session, the S4Vectors syntax was not working.
#This is when I discovered that when I wrote the script, S4Vectors' rename function had been overriding the dplyr rename function,
#and I needed to load S4Vectors to run the script with the syntax I had written.

#import normalized gene expression counts
normcounts <- read_xlsx("D_SourceDataFigs.xlsx", sheet = "Fig. 2b")

#Clean data
normcounts<- normcounts%>% #Rename plant parts to clean up names that will appear on figure
  S4Vectors::rename(
    "Fourth_Leaf" = "Fourth Leaf",
    "Sixth_Leaf" = "Sixth Leaf",
    "MatureLeaf" = "Mature Leaf",
    "EVM" = "Meristem: Veg. early",
    "MVM" = "Meristem: Veg. mid",
    "LVM" = "Meristem: Veg. late",
    "TM" = "Meristem: Transition",
    "FM" = "Meristem: Floral",
    "SYM" = "Meristem: Sympodial",
    "SIM" = "Meristem: Inflorescence",
    "FlowerBud" = "Bud",
    "OnecmFruit" = "Fruit: 1 cm",
    "GreenFruit" = "Fruit: Green",
    "BreakerFruit" = "Fruit: Breaker",
    "MatureFruit" = "Fruit: Mature"
  )%>%
  select(-c("EVM_Leaf", "MVM_Leaf", "LVM_Leaf", "TM_Leaf", "SYM_Leaf"))%>%#Remove data not used in the publication figure
  pivot_longer(#pivot data longer so it is tidy
    cols = c("Root",
             "Fourth Leaf",
             "Sixth Leaf",
             "Mature Leaf",
             "Meristem: Veg. early",
             "Meristem: Veg. mid",
             "Meristem: Veg. late",
             "Meristem: Transition",
             "Meristem: Inflorescence",
             "Meristem: Floral",
             "Meristem: Sympodial",
             "Bud",
             "Flower",
             "Pollen",
             "Fruit: 1 cm",
             "Fruit: Green",
             "Fruit: Breaker",
             "Fruit: Mature"
),
    names_to = "Plant_Part",
    values_to = "Norm_Counts_TPM")%>%
  mutate(Plant_Part = as.factor(Plant_Part))%>%#Convert plant part column to factor data and then relevel factors to the order they need to appear in plot
  mutate(Plant_Part = fct_relevel(Plant_Part, "Root",
                                  "Fourth Leaf",
                                  "Sixth Leaf",
                                  "Mature Leaf",
                                  "Meristem: Veg. early",
                                  "Meristem: Veg. mid",
                                  "Meristem: Veg. late",
                                  "Meristem: Transition",
                                  "Meristem: Floral",
                                  "Meristem: Inflorescence",
                                  "Meristem: Sympodial",
                                  "Bud",
                                  "Flower",
                                  "Pollen",
                                  "Fruit: 1 cm",
                                  "Fruit: Green",
                                  "Fruit: Breaker",
                                  "Fruit: Mature"
))%>%
  mutate(gene = as.factor(gene))%>%#Convert gene column to factor data and then relevel factors to the order they need to appear in plot
  mutate(gene = fct_relevel(gene, "Solyc02g083520.2.1", "Solyc02g061990.2.1"))

#Plot
normcounts%>%
  ggplot(aes(x = Plant_Part, y = Norm_Counts_TPM, fill = gene))+
  geom_col()+
  facet_wrap(vars(gene), ncol = 1, scales = "free")+#use facet_wrap() to show plots for both genes in one panel;
  #ncol = 1 makes the facet into 1 column; scales = free allows for the use of facetted_pos_scales to manually set the axis limits differently for each plot
  facetted_pos_scales(y = list(gene == "Solyc02g083520.2.1" ~ scale_y_continuous(limits = c(0, 1200)),
                                gene == "Solyc02g061990.2.1" ~ scale_y_continuous(limits = c(0, 550))))+
  theme_few()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2), legend.position = "none")+#angle, hjust, and vjust allow for customization of x axis labels positioning
  labs(x = "Plant Part", y = "Normalized counts (TPM)")+
  scale_fill_manual(values=c("gold1", "dodgerblue4"))


