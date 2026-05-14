########
setwd("~/CM_515/FinalProj/data")

#install.packages("readxl")
library(readxl)
library(tidyverse)
library(S4Vectors)
library(ggthemes)
library(dplyr)
library(tidyr)

#Import data
allfreq <- read_xlsx("D_SourceDataFigs.xlsx", sheet = "Fig. 2d")

#Clean data
allfreq2 <- allfreq %>%
  rename(
    "accession" = "Accession",
    "GT_SL2.50ch02_335611" = "Alleles",
    "group" = "Group",
    "excluded" = "Excluded"
  )%>%
  filter(is.na(Excluded))%>%#Filter out data listed as excluded
  filter(Alleles != "./.")%>%#Filter out data without information in the allele column
  filter(!Group %in% c("Mixture", "SLC_SP_Peru", "SLC_non_Andean", "unk"))%>% #Filter out additional data excluded by authors
  separate_wider_delim(Alleles, delim = "/", names = c("Allele_1", "Allele_2"))%>%#Separate zygosity information into two columns, 1 for each allele
  pivot_longer(
    cols = c(Allele_1, Allele_2),
    names_to = "Allele_#",
    values_to = "Allele_Version"
  )%>%#With this pivot longer, each accession now has TWO rows, one for each of its two alleles (whether homozygous or heterozygous)
  mutate(Group = gsub("wild", "Wild", Group))%>%#Rename groups to clean names that will appear on plot
  mutate(Group = gsub("SP_Ecuador", "S. pimpinellifolium", Group))%>%
  mutate(Group = gsub("SP_Peru", "S. pimpinellifolium", Group))%>%
  mutate(Group = gsub("SLC_Ecuador", "S. lycopersicum var. cerasiforme Ecuador", Group))%>%
  mutate(Group = gsub("SLC_Peru", "S. lycopersicum var. cerasiforme Peru", Group))%>%
  mutate(Group = gsub("SLL_Vintage_1", "S. lycopersicum vintage", Group))%>%
  mutate(Group = gsub("SLL_Vintage_2", "S. lycopersicum vintage", Group))%>%
  mutate(Group = gsub("SLL_Vintage_3", "S. lycopersicum vintage", Group))%>%
  mutate(Group = gsub("SLL_Fresh_Vintage", "S. lycopersicum fresh", Group))%>%#must put this above SLL_Fresh rename for both to be changed
  mutate(Group = gsub("SLL_Fresh", "S. lycopersicum fresh", Group))%>%
  mutate(Group = gsub("SLL_Processing", "S. lycopersicum processing", Group))
  
########################



#Summarise the count of each allele version per group into a new column

sum_allfreq <- allfreq2 %>%
  group_by(Group, Allele_Version) %>%
  dplyr::summarise(count = n(), .groups = "drop")%>%
  pivot_wider(
    names_from = Allele_Version,
    values_from = count
  ) %>%
   rename(
     "C" = "No_of_Cs",
     "T" = "No_of_Ts"
   )

#Pivot longer so data is tidy
sum_allfreq <- sum_allfreq %>%
  pivot_longer(cols = c(No_of_Cs, No_of_Ts),
               names_to = "Allele---Version",
               values_to = "Count"
               )

#Convert NAs (many currently in Count column) to 0s
sum_allfreq[is.na(sum_allfreq)] <- 0

#Currently multiple rows per Group + Allele Version (because there were many accessions in each group).
#Summarise to combine these rows
sum_allfreq<- sum_allfreq%>%
  group_by(Group, `Allele---Version`)%>%
  dplyr::summarise(Allele_Count = sum(Count))

#Change Allele Version content to actual names of variants
#Convert Allele Version column to factor data and relevel factors to the order that they appear on plot
sum_allfreq <- sum_allfreq%>%
  mutate(`Allele---Version` = gsub("No_of_Cs", "SSP2 S169", `Allele---Version`))%>%
  mutate(`Allele---Version` = gsub("No_of_Ts", "SSP2 F169", `Allele---Version`))%>%
  mutate(`Allele---Version`= as.factor(`Allele---Version`))%>%
  mutate(Group = as.factor(Group))%>%
  mutate(`Allele---Version` = fct_relevel(`Allele---Version`, "SSP2 S169", "SSP2 F169"))%>%
  mutate(Group = fct_relevel(Group, "Wild", "Galapagos", "S. pimpinellifolium", "S. lycopersicum var. cerasiforme Ecuador", "S. lycopersicum var. cerasiforme Peru", "S. lycopersicum vintage", "S. lycopersicum fresh", "S. lycopersicum processing"
))

#make df with a column containing the number of accessions per Group (this is the n value on plot)  
accession_count_pergroup <- sum_allfreq%>%
  group_by(Group)%>%
  summarise(accession_count_pergroup = ((sum(Allele_Count))/2))

#join dataframes
joined <- full_join(sum_allfreq, accession_count_pergroup)

#combine the following columns (Group and Accession Count per Group) into a new column for labeling purposes
joined$Group_with_acc_n <- paste(joined$Group, paste0("n=", joined$accession_count_pergroup, ""), sep = "\n")


#Plot

joined%>%
  ggplot(aes(x = `Group`, y = Allele_Count, fill = `Allele---Version`))+
  geom_bar(position="fill", stat="identity")+#position = fill stacks the two bars for each group on top of each other and sets the values to % rather than absolute values;
  #stat = identity sets the height of the bars to the values assigned to y in prior ggplot aesthetics mapping
  theme_few()+
  scale_fill_manual(values = c("dodgerblue4", "seagreen3"))+#assign color to allele versions
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 10), legend.title = element_blank(), legend.position = "top")+#customize positioning of x axis labels
  labs(x = "", y = "Allele frequency (%)")+
  scale_x_discrete(labels = setNames(joined$Group_with_acc_n, joined$Group))#use scale_x_discrete to keep the x axis as Group but LABEL the x axis with the combined column Group_with_acc_n

