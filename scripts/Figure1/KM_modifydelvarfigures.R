#Figure Clean Up: Delvar Data Analysis
#Must run delvar_data-analysis.R FIRST and keep all objects, etc. loaded.


library(tidyverse)

### FIGURE 1A:
#Barplots DELETERIOUS variants

#Clean Data
#Rename column contents to wording that should appear on figure 
altcount_list$GT[altcount_list$GT == "0/1"] <- "Heterozygous"
altcount_list$GT[altcount_list$GT == "1/1"] <- "Homozygous"
altcount_list$SPECIES[altcount_list$SPECIES == "SLC"] <- "S. lycopersicum var. cerasiforme"
altcount_list$SPECIES[altcount_list$SPECIES == "SLL"] <- "S. lycopersicum var. lycopersicum"
altcount_list$SPECIES[altcount_list$SPECIES == "SP"] <- "S. pimpinellifolium"

#Remove unnecessary text from column contents
altcount_list <- altcount_list %>%
  mutate(SAMPLE = gsub("Porter_", "", SAMPLE))%>%
  mutate(SAMPLE = gsub("_.*", "", SAMPLE))

#Plot

#Modified code to plot the data as it appears in publication, with the color of bars representing zygosity,
#and a row of colored squares along the bottom representing species.

#Use two geom_bar geometries to plot (1) the bars representing the number of deleterious variants and
#(2) the row of colored squares at the bottom of the plot indicating species.
#Use scale_fill_manual to assign colors to both bar geometries.
#Modify position of legend, remove legend title.

figure_a <- altcount_list %>%
  ggplot(aes(x= reorder(SAMPLE, count))) + 
  geom_bar(aes(y= count, fill=GT), position="stack", stat="identity") +
  geom_bar(aes(y=-300, fill=SPECIES), position="stack", stat="identity") +
  theme_light(base_size = 6) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", legend.title = element_blank()) +
  scale_fill_manual(values=c("S. lycopersicum var. cerasiforme" = '#d95f02',"S. lycopersicum var. lycopersicum" = '#7570b3',"S. pimpinellifolium" = '#1b9e77', "Heterozygous" = '#9FAFC3',"Homozygous" = '#3F5F88'))+
  labs(x = "", y = "Number of variants")+
  coord_cartesian(expand = FALSE)

figure_a

ggsave(paste0(date, "_", experiment, "/barplot_ALTGT_counts_per_accession_GTANDSPECIES.pdf"), width = 5, height = 5, useDingbats=FALSE)

####################################################################


# FIGURE 1 B and D:
# plot SIFT_SCORE

#Made minor modifications to code to clean up appearance of plots and adjust them to be similar to publication figures,
#including adjusting location of legend, adding labels, setting theme to classic, and using ggarrange()
#to arrange the plots into one panel

p1 <-  ggplot(CETS_dat.var.counts) +
  geom_point(aes(reorder(GENE_ID, count, mean), y=SIFT_SCORE, size = count, color = -SIFT_SCORE)) +
  theme_light() +
  scale_y_reverse() +
  coord_flip() +
  theme(axis.title.y = element_blank(), axis.text.x = element_text(angle = 0)) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", linewidth=1) +#Changed size to linewidth because of following warning "Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0. Please use `linewidth` instead."
  theme(legend.position = "right", legend.justification = c(0, 1)) +
  scale_colour_viridis() +
  ggtitle("CETS")+
  labs(x = "", y = "SIFT Score", size = "Count", color = "SIFT Score")+
  theme_classic()

  
#p1

p2 <- ggplot(bZIPA_dat.var.counts) +
  geom_point(aes(reorder(GENE_ID, count, mean), y=SIFT_SCORE, size = count, color = -SIFT_SCORE)) +
  theme_light() +
  scale_y_reverse() +
  coord_flip() +
  theme(axis.title.y = element_blank(), axis.text.x = element_text(angle = 0)) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", linewidth=1) +#Changed size to linewidth because of following warning "Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0. Please use `linewidth` instead."
  theme(legend.position = "right") +
  scale_colour_viridis() +
  ggtitle("bZIP")+
  labs(x = "", y = "SIFT Score", size = "Count", color = "SIFT Score")+
  theme_classic()

  
#p2
figure_bd <- ggarrange(p1, p2,
                    labels = c("b", "d"),
                    ncol = 1, nrow = 3,
                    legend = "right",
                    common.legend = TRUE,
                    heights = c(1, 1.5))
figure_bd



#######################################


#FIGURE 1 C AND E


#Made minor modifications to code to clean up appearance of plots and adjust them to be similar to publication figures,
#including adding titles, adding labels, setting the theme to classic, and using ggarrange() to arrange both plots
#into one panel.

# deleterious only
CETS_counts_del <- dat_CETS %>%
  filter(!grepl("TOLERATED", SIFT_PREDICTION)) %>%
  group_by(NONSYN_ID, SPECIES, GT) %>% 
  dplyr::summarise(COUNT = n())

CE <- CETS_counts_del %>%
  filter(GT %in% c("1/1", "0/1")) %>%
  ggplot(aes(fill=SPECIES, y=COUNT, x = NONSYN_ID)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~factor(SPECIES, levels=c('SP','SLC','SLL')), ncol=3) +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top") +
  theme_light(base_size = 6) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c('#d95f02','#7570b3','#1b9e77'))+
  ggtitle("CETS")+
  labs(y = "Number of accessions", x = "")+
  theme_classic()


# all deleterious
bZIPA_counts_del <- dat_bZIPA %>%
  filter(!grepl("TOLERATED", SIFT_PREDICTION)) %>%
  group_by(NONSYN_ID, SPECIES, GT) %>% 
  dplyr::summarise(COUNT = n())

bZ <- bZIPA_counts_del %>%
  filter(GT %in% c("1/1", "0/1")) %>%
  ggplot(aes(fill=SPECIES, y=COUNT, x = NONSYN_ID)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~factor(SPECIES, levels=c('SP','SLC','SLL')), ncol=3) +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top") +
  theme_light(base_size = 6) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c('#d95f02','#7570b3','#1b9e77'))+
  ggtitle("bZIP")+
  labs(y = "Number of accessions", x = "")+
  theme_classic()

figure_ce <- ggarrange(CE, bZ,
                    labels = c("c", "e"),
                    ncol = 1, nrow = 3,
                    legend = "none",
                    heights = c(1, 1.5))

figure_ce