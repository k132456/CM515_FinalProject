#!/bin/bash
## =============================================================================
## fig3C_motif_analysis.sh
##
## Produces Figure 3C from Glaus et al. 2025 Nat Genet
## De novo motif logos for SSP2-S169 (SpimSSP2) and SSP2-F169 (SlycSSP2)
## private peaks using MEME Suite.
##
## INPUTS REQUIRED:
##   - SpimSSP2_peaks_private_fc3_200bp.bed   (from dapseq-analysis.R)
##   - SlycSSP2_peaks_private_fc3_200bp.bed   (from dapseq-analysis.R)
##   - Sweet-100 v2.0 genome FASTA
##
## DEPENDENCIES (install once, see Step 1 below):
##   conda install -c bioconda bedtools meme -y
##
## USAGE:
##   1. Open Ubuntu terminal
##   2. conda activate deeptools
##   3. bash /mnt/c/Users/kmark/Documents/CM_515/fig3C_motif_analysis.sh
## =============================================================================


## =============================================================================
## 0. SETUP — adjust these paths to match your file locations
## =============================================================================

## Input BED files (from dapseq-analysis.R output folder)
## Use the Windows path converted to WSL format (/mnt/c/...)
BED_SPIM="/mnt/c/Users/kmark/Documents/CM_515/SpimSSP2_peaks_private_fc3_200bp.bed"
BED_SLYC="/mnt/c/Users/kmark/Documents/CM_515/SlycSSP2_peaks_private_fc3_200bp.bed"

## Genome FASTA
GENOME="/mnt/c/Users/kmark/Documents/CM_515/SollycSweet-100_v2.0.fasta"

## Output directory
OUT_DIR="/mnt/c/Users/kmark/Documents/CM_515/fig3C_output"
mkdir -p "${OUT_DIR}"


## =============================================================================
## STEP 1: INSTALL DEPENDENCIES (run once)
## If bedtools and meme are already installed, skip this step.
## =============================================================================

echo "============================================="
echo "Step 1: Checking dependencies..."
echo "============================================="

if ! command -v bedtools &> /dev/null; then
    echo "  bedtools not found — installing..."
    conda install -c bioconda bedtools -y
else
    echo "  bedtools found: $(bedtools --version)"
fi

if ! command -v meme &> /dev/null; then
    echo "  MEME not found — installing..."
    conda install -c bioconda meme -y
else
    echo "  MEME found: $(meme --version)"
fi


## =============================================================================
## STEP 2: FIX WINDOWS LINE ENDINGS IN BED FILES
## BED files created on Windows may have \r\n line endings that break
## bedtools on Linux. This step removes them.
## =============================================================================

echo ""
echo "============================================="
echo "Step 2: Fixing line endings in BED files..."
echo "============================================="

sed -i 's/\r//' "${BED_SPIM}"
sed -i 's/\r//' "${BED_SLYC}"

echo "  SpimSSP2 BED: $(wc -l < ${BED_SPIM}) peaks"
echo "  SlycSSP2 BED: $(wc -l < ${BED_SLYC}) peaks"
echo "  Preview of SpimSSP2 BED:"
head -3 "${BED_SPIM}"


## =============================================================================
## STEP 3: EXTRACT DNA SEQUENCES USING BEDTOOLS GETFASTA
## For each BED file, extract the corresponding DNA sequences from the
## genome FASTA. The -s flag respects strand orientation.
## =============================================================================

echo ""
echo "============================================="
echo "Step 3: Extracting peak sequences with bedtools..."
echo "============================================="

FASTA_SPIM="${OUT_DIR}/SpimSSP2_private_sequences.fasta"
FASTA_SLYC="${OUT_DIR}/SlycSSP2_private_sequences.fasta"

bedtools getfasta \
    -fi "${GENOME}" \
    -bed "${BED_SPIM}" \
    -fo "${FASTA_SPIM}"
echo "  SpimSSP2 sequences extracted: $(grep -c '>' ${FASTA_SPIM})"

bedtools getfasta \
    -fi "${GENOME}" \
    -bed "${BED_SLYC}" \
    -fo "${FASTA_SLYC}"
echo "  SlycSSP2 sequences extracted: $(grep -c '>' ${FASTA_SLYC})"


## =============================================================================
## STEP 4: RUN MEME FOR DE NOVO MOTIF DISCOVERY
## Parameters match exactly those stated in the paper Methods section:
##   -dna           : input sequences are DNA
##   -mod zoops     : zero or one occurrence per sequence model
##   -nmotifs 3     : find top 3 motifs
##   -minw 6        : minimum motif width 6bp
##   -maxw 15       : maximum motif width 15bp
##   -maxsites 1000 : maximum number of sites
##   -objfun classic: classic objective function
##   -revcomp       : search both strands
##   -markov_order 0: 0th order background model
##
## Note: SpimSSP2 uses top 1000 peaks (already selected in dapseq-analysis.R)
##       SlycSSP2 uses all 426 private peaks (fewer than 1000 so use all)
## =============================================================================

echo ""
echo "============================================="
echo "Step 4: Running MEME for SpimSSP2 (SSP2-S169)..."
echo "============================================="

MEME_OUT_SPIM="${OUT_DIR}/meme_SpimSSP2"

meme "${FASTA_SPIM}" \
    -dna \
    -mod zoops \
    -nmotifs 3 \
    -minw 6 \
    -maxw 15 \
    -maxsites 1000 \
    -objfun classic \
    -revcomp \
    -markov_order 0 \
    -oc "${MEME_OUT_SPIM}" \
    -p 4

echo "  MEME complete for SpimSSP2."
echo "  Results in: ${MEME_OUT_SPIM}"

echo ""
echo "============================================="
echo "Step 4b: Running MEME for SlycSSP2 (SSP2-F169)..."
echo "============================================="

MEME_OUT_SLYC="${OUT_DIR}/meme_SlycSSP2"

meme "${FASTA_SLYC}" \
    -dna \
    -mod zoops \
    -nmotifs 3 \
    -minw 6 \
    -maxw 15 \
    -maxsites 1000 \
    -objfun classic \
    -revcomp \
    -markov_order 0 \
    -oc "${MEME_OUT_SLYC}" \
    -p 4

echo "  MEME complete for SlycSSP2."
echo "  Results in: ${MEME_OUT_SLYC}"


## =============================================================================
## STEP 5: SUMMARISE OUTPUTS
## MEME produces several output files. The most relevant for Figure 3C are:
##   meme.html     : interactive results viewer (open in browser)
##   meme.txt      : plain text results
##   logo1.png     : sequence logo for top motif  <-- THIS IS FIGURE 3C
##   logo2.png     : sequence logo for 2nd motif
##   logo3.png     : sequence logo for 3rd motif
## =============================================================================

echo ""
echo "============================================="
echo "Step 5: Output summary"
echo "============================================="

echo ""
echo "SpimSSP2 (SSP2-S169) output files:"
ls -lh "${MEME_OUT_SPIM}/" | grep -v "^total" | awk '{print "  "$NF"  ("$5")"}'

echo ""
echo "SlycSSP2 (SSP2-F169) output files:"
ls -lh "${MEME_OUT_SLYC}/" | grep -v "^total" | awk '{print "  "$NF"  ("$5")"}'

echo ""
echo "============================================="
echo "Figure 3C logos are:"
echo "  ${MEME_OUT_SPIM}/logo1.png  (SSP2-S169 top motif)"
echo "  ${MEME_OUT_SLYC}/logo1.png  (SSP2-F169 top motif)"
echo ""
echo "Open meme.html in your browser for the full"
echo "interactive results including all motifs."
echo "============================================="
echo ""
echo "Done!"
