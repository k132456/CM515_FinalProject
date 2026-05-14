#This script contains a modified version of Section 2 of dapseq-analysis.R (meant to replace the original Section 2), with modifications made in order to 
#clean up figure appearance for Figures 3A and 3D. These modifications exist between lines 85 and 126 and between lines 195 and 221 of THIS script.

#85-126: modifications clean up labels and title for 3D
#195-221: modifications clean up labels and title for 3A


# 2. VISUALISATION OF SIGNIFICANT PEAKS AND ASSOCIATED GENES

# 2.1 Compare gene and peak sets in Venn-Euler diagrams to determine the intersection of peaks/genes bound by the different bait proteins. Genes of different sectors are written with gene descriptions into .csv files.

# 2.1.1 Gene comparisons
type="genes"
# Set data
SSP <- SSP_genes
SlycSSP2 <- SlycSSP2_genes
SpimSSP2 <- SpimSSP2_genes

# Calculate total number of unique features
length(SSP)
length(SlycSSP2)
length(SpimSSP2)
combined <- c(SSP,SlycSSP2,SpimSSP2)
length(combined)
length(unique(combined))
total <- length(unique(combined))

# Calculate overlap between vectors
overlap <- calculate.overlap(
  x = list(
    "SSP" = SSP,
    "SlycSSP2" = SlycSSP2,
    "SpimSSP2" = SpimSSP2
  )
);

# rename the overlap lists (for a 3-way overlap)
names(overlap) <- c("a123", "a12", "a13", "a23", "a1", "a2", "a3")

# calculate length of overlaps
a1=length(overlap$a1)
a2=length(overlap$a2)
a3=length(overlap$a3)
a12=length(overlap$a12)
a13=length(overlap$a13)
a23=length(overlap$a23)
a123=length(overlap$a123)

# write overlaps with gene annotation to file
overlap_a1 <- as_data_frame(overlap$a1) %>%
  rename(geneid=value) %>%
  left_join(desc, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SSP_fc",fc,".csv")))

overlap_a2 <- as_data_frame(overlap$a2) %>%
  rename(geneid=value) %>%
  left_join(desc, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SlycSSP2_fc",fc,".csv")))

overlap_a3 <- as_data_frame(overlap$a3) %>%
  rename(geneid=value) %>%
  left_join(desc, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SpimSSP2_fc",fc,".csv")))

overlap_a12 <- as_data_frame(overlap$a12) %>%
  rename(geneid=value) %>%
  left_join(desc, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SSP&SlycSSP2_fc",fc,".csv")))

overlap_a13 <- as_data_frame(overlap$a13) %>%
  rename(geneid=value) %>%
  left_join(desc, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SSP&SpimSSP2_fc",fc,".csv")))

overlap_a23 <- as_data_frame(overlap$a23) %>%
  rename(geneid=value) %>%
  left_join(desc, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SlycSSP2&SpimSSP2_fc",fc,".csv")))

overlap_a123 <- as_data_frame(overlap$a123) %>%
  rename(geneid=value) %>%
  left_join(desc, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SSP&SlycSSP2&SpimSSP2_fc",fc,".csv")))

# Make Euler object
euler <- euler(c("SSP" = a1,#Rename proteins to clean names that should appear in file
                 "SSP2F169" = a2,
                 "SSP2S169" = a3,
                 "SSP&SSP2F169" = a12,
                 "SSP&SSP2S169" = a13,
                 "SSP2F169&SSP2S169" = a23,
                 "SSP&SSP2F169&SSP2S169" = a123))

# Select colors 
myCol <- c("#FCC427", "#35B779", "#31688E") # Custom for SSP-SlycSSP2-SpimSSP2

length(overlap$a1)
SSP_genes_unique <- overlap$a1
length(overlap$a2)
SlycSSP2_genes_unique <- overlap$a2
length(overlap$a3)
SpimSSP2_genes_unique <- overlap$a3
length(overlap$a13)
SSP_SpimSSP2_genes_shared <- overlap$a13
length(overlap$a13)
SSP_SpimSSP2_genes_shared <- overlap$a13
length(overlap$a123)
SSP_SlSpSSP2_genes_shared <- overlap$a123


#######################
# Plot Euler object
library(grid)
euler.plot <- plot(euler,
                   quantities = list(type = c("counts")),
                   fills = myCol,
                   alpha=0.5, # change transparency depending on the used color palette
                   edges=F,
                   main="Genes\nTotal = 8,395",#Clean up title
                   labels = FALSE)#Remove labels and instead plot with grid.text (below)
euler.plot

grid.text("SSP", x = 0.16, y = 0.78, gp = gpar(fontsize = 12))#Label circles with grid.text instead of within plot() in order to control label location (with x= and y=) and font size.
grid.text("SSP2S169", x = 0.81, y = 0.81, gp = gpar(fontsize = 12))
grid.text("SSP2F169", x = 0.93, y = 0.26, gp = gpar(fontsize = 12))


# Save plot to file
# 1.Open a pdf file
pdf((file.path(paste0(date, "_", experiment, "/Venn/Venn_", type, "_fc", fc, ".pdf"))), width = 5, height = 4, useDingbats=FALSE)
# 2. Create a plot
euler.plot
# Close the pdf file
dev.off() 


# 2.1.2 Peak comparisons
type="peaks"
# Set data
SSP <- SSP_peaks
SlycSSP2 <- SlycSSP2_peaks
SpimSSP2 <- SpimSSP2_peaks

# Calculate total number of unique features
length(SSP)
length(SlycSSP2)
length(SpimSSP2)
combined <- c(SSP,SlycSSP2,SpimSSP2)
length(combined)
length(unique(combined))
total <- length(unique(combined))
# Calculate overlap between vectors
overlap <- calculate.overlap(
  x = list(
    "SSP" = SSP,
    "SlycSSP2" = SlycSSP2,
    "SpimSSP2" = SpimSSP2
  )
);

# rename the overlap lists (for a 3-way overlap)
names(overlap) <- c("a123", "a12", "a13", "a23", "a1", "a2", "a3")
# calculate length of overlaps
a1=length(overlap$a1)
a2=length(overlap$a2)
a3=length(overlap$a3)
a12=length(overlap$a12)
a13=length(overlap$a13)
a23=length(overlap$a23)
a123=length(overlap$a123)

# write overlaps with gene annotation to file
overlap_a1 <- as_data_frame(overlap$a1) %>%
  rename(peak_id=value) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SSP_fc",fc,".csv")))
overlap_a2 <- as_data_frame(overlap$a2) %>%
  rename(peak_id=value) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SlycSSP2_fc",fc,".csv")))
overlap_a3 <- as_data_frame(overlap$a3) %>%
  rename(peak_id=value) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SpimSSP2_fc",fc,".csv")))
overlap_a12 <- as_data_frame(overlap$a12) %>%
  rename(peak_id=value) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SSP&SlycSSP2_fc",fc,".csv")))
overlap_a13 <- as_data_frame(overlap$a13) %>%
  rename(peak_id=value) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SSP&SpimSSP2_fc",fc,".csv")))
overlap_a23 <- as_data_frame(overlap$a23) %>%
  rename(peak_id=value) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SlycSSP2&SpimSSP2_fc",fc,".csv")))
overlap_a123 <- as_data_frame(overlap$a123) %>%
  rename(peak_id=value) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SSP&SlycSSP2&SpimSSP2_fc",fc,".csv")))

# Make Euler object
euler <- euler(c("SSP" = a1,#Rename proteins to clean names that should appear in file
                 "SSP2F169" = a2,
                 "SSP2S169" = a3,
                 "SSP&SSP2F169" = a12,
                 "SSP&SSP2S169" = a13,
                 "SSP2F169&SSP2S169" = a23,
                 "SSP&SSP2F169&SSP2S169" = a123))

# Select colors 
myCol <- c("#FCC427", "#35B779", "#31688E") # Custom for SSP-SlycSSP2-SpimSSP2

# Plot Euler object
library(grid)
euler.plot <- plot(euler,
                   quantities = list(type = c("counts")),
                   fills = myCol,
                   alpha=0.5, # change transparency depending on the used color palette
                   edges=F,
                   main="Peaks\nTotal = 14,091",#Clean up title
                   labels = FALSE)#Remove labels and instead plot with grid.text (below)
euler.plot

grid.text("SSP", x = 0.16, y = 0.78, gp = gpar(fontsize = 12))#Label circles with grid.text instead of within plot() in order to control label location (with x= and y=) and font size.
grid.text("SSP2S169", x = 0.81, y = 0.81, gp = gpar(fontsize = 12))
grid.text("SSP2F169", x = 0.93, y = 0.245, gp = gpar(fontsize = 12))



######################


# Save plot to file
# 1.Open a pdf file
pdf((file.path(paste0(date, "_", experiment, "/Venn/Venn_", type, "_fc", fc, ".pdf"))), width = 5, height = 4, useDingbats=FALSE)
# 2. Create a plot
euler.plot
# Close the pdf file
dev.off()
