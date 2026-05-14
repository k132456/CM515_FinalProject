# AnnotationForge
# build custom gene ontology databases
# https://bioconductor.org/packages/release/bioc/html/AnnotationForge.html


#The authors imported data in the original script via a file ("ITAG4.0.descriptions.shortid"),
#that they did not provide and was not available on the ITAG4.0 repository (neither under the same or a different file name).
#Based on the dataframes they made from the information in this file, I was able to determine what data was used
#and import a different file from the ITAG4.0 repository that contained this data.
#I called the package rtracklayer in order to import the data from this file, which contained gene model data in
#the form of a GFF file.
#I created a dataframe (gffdat_genes_df) that contained the critical information contained in the corresponding
#dataframe that the authors had created from data in "ITAG4.0.descriptions.shortid".
#(Note: In creating this dataframe, I completed some tasks (like renaming columns) that in the original script were done in
#subsequent subsetted dataframes, so I removed these redundant tasks from the later code chunks.)
#I was able to then run the rest of the script as is.

setwd("~/CM_515/FinalProj/data")
home <- getwd()
#BiocManager::install("AnnotationForge")

# Load libraries
library(AnnotationForge)
library(tidyverse)
library(rtracklayer)

##################################################################################################
################Code in the section below has been added or changed###############################
#Import gff

gffdat<- import("ITAG4.0_gene_models.gff")
 
gffdat_genes <- gffdat[gffdat$type == "gene"]
 
gffdat_genesdf<- as.data.frame(gffdat_genes)

#gff edits
gffdat_genesdf <- gffdat_genesdf[, -c(2:9, 11, 13:21)]
gffdat_genesdf$ID <- gsub("gene:", "", gffdat_genesdf$ID)
colnames(gffdat_genesdf) <- c("CHROMOSOME", "GID", "GENENAME")
gffdat_genesdf$CHROMOSOME <- gsub("SL4.0", "", gffdat_genesdf$CHROMOSOME)
gffdat_genesdf <- gffdat_genesdf%>%
  relocate("CHROMOSOME", .after = "GENENAME")


# Make Sym df
#Changes made here are a result of importing information directly into a GID column in gffdat_genesdf
#rather than having both a SYMBOL and a GID column as the original script did.
Sym <- gffdat_genesdf[,c(1,2)]
head(Sym)
colnames(Sym) <- c("GID", "GENENAME")
Sym$GID <- gsub("\\..*","",Sym$GID)
length(unique(Sym$GID))
#length(unique(Sym$SYMBOL))
#Sym$GID <- Sym$SYMBOL
head(Sym)
tail(Sym)
nrow(Sym)
length(unique(Sym$GID))

# Make Chr df
#head(descdat)
Chr <- gffdat_genesdf[,c(1,3)]
head(Chr)
#colnames(Chr) <- c("GID", "CHROMOSOME") #Did this above
#Chr$CHROMOSOME <- gsub("SL4.0","",Chr$CHROMOSOME) #Did this above

length(unique(Chr$GID))
Chr$GID <- gsub("\\..*","",Chr$GID)
length(unique(Chr$GID))

head(Chr)
tail(Chr)

################Code in the section above has been added or changed###############################
##################################################################################################

godat <- read_delim("go.sly.csv", skip = 8, delim="\t") # downloaded from https://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_05/GO/go.sly.csv.gz
head(godat)
nrow(godat)

# Make GO df
head(godat)
sum(duplicated(godat))

GO <- godat[,c(1,3,4)]
head(GO)
sum(duplicated(GO)) # number of duplicated rows
colnames(GO) <- c("GID","GO", "EVIDENCE")
nrow(GO)
GO <- GO %>% drop_na(GO) # remove NAs
nrow(GO)

length(unique(GO$GID))
GO$GID <- gsub("\\..*","",GO$GID)
length(unique(GO$GID))

head(GO)
tail(GO)

#################
# check files dfs
nrow(Sym)
nrow(Chr)
nrow(GO)

str(Sym)
str(Chr)
str(GO)

head(Sym)
head(Chr)
head(GO)


# check overlaps between Sym and GO dataframes
SymGO_overlap <- Sym %>%
  filter(GID %in% GO$GID)
nrow(SymGO_overlap)  # number of genes WITH annotation
nrow(Sym) - nrow(SymGO_overlap) # number of genes WITHOUT annotation

Sym_noGO <- Sym %>%
  filter(!GID %in% GO$GID)
nrow(Sym_noGO)
head(Sym_noGO) # genes WITHOUT annotation
tail(Sym_noGO) # genes WITHOUT annotation


# check overlaps between Chr and GO dataframes
ChrGO_overlap <- Chr %>%
  filter(GID %in% GO$GID)
nrow(ChrGO_overlap)
nrow(Chr)
nrow(ChrGO_overlap) # number of genes WITH annotation
nrow(Chr) - nrow(ChrGO_overlap) # number of genes WITHOUT annotation

Chr_noGO <- Chr %>%
  filter(!GID %in% GO$GID)
nrow(Chr_noGO)
head(Chr_noGO) # genes WITHOUT annotation
tail(Chr_noGO) # genes WITHOUT annotation

# check overlaps between GO and Sym dataframes
GOSym_overlap <- GO %>%
  filter(GID %in% Sym$GID)
nrow(GOSym_overlap) 
nrow(GO)

noGOSym_overlap <- GO %>%
  filter(!GID %in% Sym$GID)
nrow(noGOSym_overlap) # shoud be 0
nrow(GO)

noGOChr_overlap <- GO %>%
  filter(!GID %in% Chr$GID)
nrow(noGOChr_overlap) # should be 0
nrow(GO)

sum(duplicated(GO))
sum(duplicated(Sym))
sum(duplicated(Chr))

GO <- distinct(GO) # remove duplicated rows

# double-check overlaps between Sym and GO dataframes
SymGO_overlap_distinct <- Sym %>%
  filter(GID %in% GO$GID)
nrow(SymGO_overlap)  # number of genes WITH annotation before removing duplicate rows
nrow(SymGO_overlap_distinct)  # number of genes WITH annotation after removing duplicate rows

# Stats:
nrow(Sym) # 34075 (number of all genes)
nrow(SymGO_overlap) # 25750 (number of annotated genes)
nrow(GO) # 1458692 (number of GO-terms; multiple terms per gene possible)



## Call function to generate package (change package name)
makeOrgPackage(gene_info=Sym, chromosome=Chr, go=GO,
# makeOrgPackage(gene_info=Sym_dist, chromosome=Chr_dist, go=GO_dist,
               version="0.1",
               maintainer="ssoyk <sebastian.soyk@unil.ch>",
               author="ssoyk <sebastian.soyk@unil.ch>",
               outputDir = ".",
               tax_id="4081",
               genus="Solanum",
               species="lycopersicumSL40VIBS", # DO NOT USE UNDERSCORE IN NAME
               goTable="go")

## then you can call install.packages based on the return value
install.packages("./org.SlycopersicumSL40VIBS.eg.db", repos = NULL, type = "source") # = "g(genus)species" DO NOT USE UNDERSCORE IN NAME

