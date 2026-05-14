#Written by Claude AI, Checked and Modified by KM:
## =============================================================================
## Figure_3B.R
##
## Standalone script to reproduce Figure 3B from Glaus et al. 2025 Nat Genet
##
## Figure 3B: Genomic distribution bar chart showing the percentage of
## DAP-seq peaks falling in different genomic features (promoter, exon,
## intron, etc.) for SSP, SSP2-F169 (SlycSSP2), and SSP2-S169 (SpimSSP2).
##
## INPUTS REQUIRED:
##   - peak_table_fdr_0.01.tsv        (from Zenodo 13787654)
##   - SollycSweet-100_genes_v2.1.1.gff3  (genome annotation)
##
## NO bash operations required. Runs entirely in R.
##
## DEPENDENCIES:
## BiocManager::install(c("ChIPseeker", "GenomicFeatures", "GenomicRanges"))
##   install.packages(c("ggplot2", "dplyr"))
## =============================================================================


## =============================================================================
## 0. SETUP — adjust these paths
## =============================================================================
setwd("~/CM_515/FinalProj/data")
PEAK_TABLE <- "~/CM_515/FinalProj/data/peak_table_fdr_0.01.tsv"
GFF3_FILE  <- "~/CM_515/FinalProj/data/SollycSweet-100_genes_v2.1.1.gff3"
OUT_DIR    <- "~/CM_515/FinalProj/data/Figure_3B_output"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## Thresholds matching the paper
FC_CUTOFF  <- 3      # log2 fold-change cutoff vs input
TSS_UP     <- 3000   # bp upstream of TSS for promoter definition
TSS_DN     <- 2000   # bp downstream of TSS


## =============================================================================
## 1. LOAD PACKAGES
## =============================================================================

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
})
library(AnnotationDbi)

## =============================================================================
## 2. LOAD AND CLASSIFY PEAKS
## =============================================================================

message("Loading peak table...")
pt <- read.delim(PEAK_TABLE, stringsAsFactors = FALSE)

## Fix the FDR.SPP typo present in the deposited file header
if ("FDR.SPP" %in% colnames(pt) && !"FDR.SSP" %in% colnames(pt)) {
  pt <- dplyr::rename(pt, FDR.SSP = FDR.SPP)
  message("  Renamed FDR.SPP -> FDR.SSP")
}

message(sprintf("  Total peaks in table: %d", nrow(pt)))

## Subset peaks by protein using logFC cutoff (FDR filter already applied
## at table generation — all rows have FDR <= 0.01 by construction)
peaks_SSP   <- pt %>% filter(logFC.SSP      >= FC_CUTOFF)
peaks_F169  <- pt %>% filter(logFC.SlycSSP2 >= FC_CUTOFF)  # domesticated SSP2
peaks_S169  <- pt %>% filter(logFC.SpimSSP2 >= FC_CUTOFF)  # ancestral SSP2

message(sprintf("  SSP peaks:       %d", nrow(peaks_SSP)))
message(sprintf("  SSP2-F169 peaks: %d  (SlycSSP2, domesticated)", nrow(peaks_F169)))
message(sprintf("  SSP2-S169 peaks: %d  (SpimSSP2, ancestral)",    nrow(peaks_S169)))


## =============================================================================
## 3. CONVERT TO GRANGES
## =============================================================================

to_granges <- function(df) {
  GRanges(
    seqnames = df$chr,
    ranges   = IRanges(start = df$start, end = df$end)
  )
}

gr_SSP  <- to_granges(peaks_SSP)
gr_F169 <- to_granges(peaks_F169)
gr_S169 <- to_granges(peaks_S169)


## =============================================================================
## 4. BUILD TxDb FROM GFF3
##    Cached as an RDS file so it only needs to be built once. (Changed to saveDb)
## =============================================================================

txdb_cache <- file.path(OUT_DIR, "TxDb_SollycSweet100.rds")

if (file.exists(txdb_cache)) {
  message("Loading cached TxDb...")
  txdb <- readRDS(txdb_cache)
} else {
  message("Building TxDb from GFF3 (this takes a few minutes)...")
  txdb <- txdbmaker::makeTxDbFromGFF(#Added txdbmaker:: in front of makeTxDbFromGFF() because original script from Claude was using an out of date package for this
    GFF3_FILE,
    format     = "gff3",
    dataSource = "SollycSweet-100 v2.1.1",
    organism   = "Solanum lycopersicum"
  )
  saveDb(txdb, txdb_cache)#Changed to saveDb from saveRDS because saveRDS is out-of-date syntax
  message("  TxDb built and cached.")
}


## =============================================================================
## 5. ANNOTATE PEAKS WITH CHIPSEEEKER
##    annotatePeak classifies each peak into a genomic feature category
##    (promoter, 5'UTR, exon, intron, 3'UTR, downstream, intergenic)
##    based on its proximity to the nearest gene in the TxDb.
## =============================================================================

message("Annotating peaks...")

anno_SSP  <- annotatePeak(gr_SSP,  TxDb = txdb,
                          tssRegion = c(-5000, TSS_DN), verbose = FALSE)
anno_F169 <- annotatePeak(gr_F169, TxDb = txdb,
                          tssRegion = c(-5000, TSS_DN), verbose = FALSE)
anno_S169 <- annotatePeak(gr_S169, TxDb = txdb,
                          tssRegion = c(-5000, TSS_DN), verbose = FALSE)

message("  Annotation complete.")


#Sections 6 and 7 (Replacing Claude's original 6 and 7 with this updated code):

## Colours matching the paper's feature category scheme (Fig 3b)
## These match the ChIPseeker default category order:
## Promoter, 5'UTR, 3'UTR, Exon, Intron, Downstream, Intergenic
feature_colours <- c(#Changed  colors manually because Claude's color choices were very different from publication figure.
  "Promoter (<=1kb)"  = "lightskyblue",
  "Promoter (1-2kb)"  = "royalblue",
  "Promoter (2-5kb)"  = "palegreen",
  "5' UTR"            = "firebrick2",
  "3' UTR"            = "lightgoldenrod1",
  "Exon"        = "darkorange1",
  "Intron"        = "peachpuff",
  "Downstream (<=300)"= "saddlebrown",
  "Distal intergenic" = "#999999"
)

## Extract annotation statistics for each protein
get_anno_stats <- function(anno_obj, label) {
  df <- as.data.frame(anno_obj@annoStat)
  df$protein <- label
  df
}

stats_df <- bind_rows(
  get_anno_stats(anno_SSP,  "SSP"),
  get_anno_stats(anno_F169, "SSP2-F169"),
  get_anno_stats(anno_S169, "SSP2-S169")
)

stats_df <- stats_df %>%#Group promoters from three bins into one (2-5kb).
  mutate(Feature = case_when(
    Feature %in% c("Promoter (2-3kb)", 
                   "Promoter (3-4kb)", 
                   "Promoter (4-5kb)") ~ "Promoter (2-5kb)",
    TRUE ~ Feature
  )) %>%
  group_by(protein, Feature) %>%
  summarise(Frequency = sum(Frequency), .groups = "drop")

## Set protein order to match paper (SSP top, then F169, then S169)
stats_df$protein <- factor(stats_df$protein,
                           levels = c("SSP", "SSP2-F169", "SSP2-S169"))

################################
#ADDED

#Rename features to names that should show up on plot:
stats_df<-stats_df%>%
  mutate(Feature = gsub("1st Exon", "Exon", Feature))%>%
  mutate(Feature = gsub("Other Exon", "Exon", Feature))%>%
  mutate(Feature = gsub("1st Intron", "Intron", Feature))%>%
  mutate(Feature = gsub("Other Intron", "Intron", Feature))%>%
  mutate(Feature = gsub("Distal Intergenic", "Distal intergenic", Feature))


#Set feature order to match paper 
stats_df$Feature <- factor(stats_df$Feature,
                           levels = c("Distal intergenic",
  "Downstream (<=300)",
  "3' UTR",
  "Intron",
  "Exon",
  "5' UTR",
  "Promoter (2-5kb)",
  "Promoter (1-2kb)",
  "Promoter (<=1kb)"))
  
##################################

## Plot
fig3b <- ggplot(stats_df, aes(x = protein, y = Frequency, fill = Feature)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = feature_colours, na.value = "grey80") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +#sets values represented in bar to proportions not absolute values
  coord_flip() +
  labs(
    title = "",
    x     = NULL,
    y     = "Percentage of peaks",
    fill  = "Feature"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(direction = "vertical"))

print(fig3b)
ggsave(file.path(OUT_DIR, "Fig3B_genomic_distribution.pdf"),
       fig3b, width = 7, height = 3.5)
write.csv(stats_df, file.path(OUT_DIR, "Fig3B_genomic_distribution_data.csv"),
          row.names = FALSE)