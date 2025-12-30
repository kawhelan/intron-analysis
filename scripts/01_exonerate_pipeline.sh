#!/usr/bin/env bash
# ============================================================
# Exonerate → GFF3 cleanup → intron table → peptide MSA
# ------------------------------------------------------------
# What this does:
#   1) For each locus, align a peptide to its genomic sequence with exonerate
#   2) Clean exonerate output and convert to valid GFF3 (AGAT + gffread)
#   3) Rename features to consistent IDs based on the locus name
#   4) Compute intron coordinates + phases + codon index before intron
#   5) Combine introns for all loci into one TSV
#   6) Build a peptide FASTA and run MAFFT to generate a protein alignment
#
# What you need installed/available in PATH:
#   - exonerate
#   - AGAT scripts: agat_convert_sp_gxf2gxf.pl, agat_sp_add_introns.pl
#   - gffread
#   - mafft
#
# Input file naming expected per locus (you can change these patterns):
#   <IN_DIR>/<LOCUS>_peptide.fasta
#   <IN_DIR>/<LOCUS>_genomic.fasta
#
# Output structure:
#   <OUT_BASE>/<LOCUS>/... (per-locus files)
#   <OUT_BASE>/combined_introns.tsv
#   <MSA_OUT_DIR>/peptides.aln.fasta
# ============================================================

set -euo pipefail

# ----------------------------
# 0) USER SETTINGS (EDIT THESE)
# ----------------------------

# Folder that contains your per-locus peptide/genomic FASTAs
IN_DIR="data/inputs"

# Where to write all exonerate + GFF3 + intron outputs
OUT_BASE="results/exonerate_batch"

# A whitespace-separated list of locus IDs to process
# Example: "geneA geneB geneC"
LOCI="locus1 locus2 locus3"

# Exonerate settings (tweak as needed)
EXO_MODEL="protein2genome"
EXO_BESTN=1
EXO_PERCENT=25
EXO_MININTRON=20
EXO_MAXINTRON=200000

# MAFFT settings (tweak as needed)
MSA_OUT_DIR="results/intron_summaries"
MAFFT_MODE=(--localpair --maxiterate 1000 --thread -1)

# ----------------------------
# 1) BASIC SAFETY CHECKS
# ----------------------------

: "${IN_DIR:?Need to set IN_DIR}"
: "${OUT_BASE:?Need to set OUT_BASE}"
: "${LOCI:?Need to set LOCI list}"

mkdir -p "$OUT_BASE"
mkdir -p "$MSA_OUT_DIR"

# ----------------------------
# 2) RUN EXONERATE + CLEAN GFF3
# ----------------------------
# For each locus:
#   - Run exonerate peptide vs genomic
#   - Filter the output down to target GFF records
#   - Convert to valid GFF3 and add introns
#   - Normalize feature IDs (mRNA/CDS/exon/intron) so downstream scripts are consistent

for L in $LOCI; do
  OUT="$OUT_BASE/$L"
  PEP="$IN_DIR/${L}_peptide_renamed.fasta"
  GEN="$IN_DIR/${L}_genomic_renamed.fasta"

  # Check inputs exist and are not empty
  if [ ! -s "$PEP" ]; then
    echo "ERROR: missing peptide FASTA for $L ($PEP)"
    continue
  fi
  if [ ! -s "$GEN" ]; then
    echo "ERROR: missing genomic FASTA for $L ($GEN)"
    continue
  fi

  # Make a clean per-locus output folder
  mkdir -p "$OUT"
  rm -f "$OUT"/*

  echo "Running exonerate for $L ..."

  # Run exonerate and write raw output + a log file
  exonerate \
    --model "$EXO_MODEL" \
    --query "$PEP" \
    --target "$GEN" \
    --showtargetgff yes \
    --bestn "$EXO_BESTN" \
    --percent "$EXO_PERCENT" \
    --minintron "$EXO_MININTRON" \
    --maxintron "$EXO_MAXINTRON" \
    --ryo $'RESULT\t%qi\t%ti\t%qas\t%tas\t%pi\t%ps\t%V\n' \
    > "$OUT/exo.raw.gff" \
    2> "$OUT/exo.log"

  # Remove comment lines and keep only rows that look like real GFF records
  grep -v '^#' "$OUT/exo.raw.gff" \
    | awk -F'\t' 'NF>=8 && $1!="RESULT" && $4~/^[0-9]+$/ && $5~/^[0-9]+$/' \
    > "$OUT/exo.target.gff"

  # Ensure there are 9 GFF columns (attributes column exists)
  awk -F'\t' 'BEGIN{OFS="\t"} ($3 ~ /^(gene|mRNA|exon|CDS|cds|intron)$/){
      if(NF<9)$9=".";
      print
    }' "$OUT/exo.target.gff" > "$OUT/exo.target.gff9"

  # Convert to proper GFF3
  agat_convert_sp_gxf2gxf.pl --gxf "$OUT/exo.target.gff9" -o "$OUT/exo.conv.gff3"

  # Add intron features (based on exon/CDS structure)
  agat_sp_add_introns.pl --gff "$OUT/exo.conv.gff3" -o "$OUT/exo.with_introns.gff3"

  # gffread sanity/fix pass (forces consistent formatting, checks coordinates vs genome)
  gffread -E -F -g "$GEN" -o "$OUT/exo.fixed.gff3" "$OUT/exo.with_introns.gff3"

  if [ ! -s "$OUT/exo.fixed.gff3" ]; then
    echo "ERROR: no fixed GFF3 produced for $L"
    continue
  fi

  # Grab the first mRNA ID so we can rename it consistently
  MRNA_ID=$(
    awk -F'\t' 'tolower($3)=="mrna"{
      if(match($9,/ID=([^;]+)/,m)) print m[1]
    }' "$OUT/exo.fixed.gff3" | head -n1
  )

  if [ -z "${MRNA_ID:-}" ]; then
    echo "ERROR: no mRNA feature found for $L"
    continue
  fi

  # Define new IDs we want to use
  NEW_MRNA_ID="${L}.mRNA"
  NEW_CDS_ID="${L}.cds"

  # Rewrite IDs/Parents for mRNA/CDS/exon/intron so everything points to <L>.mRNA
  awk -F'\t' -v OFS='\t' \
    -v OLD="$MRNA_ID" \
    -v NM="$NEW_MRNA_ID" \
    -v NC="$NEW_CDS_ID" \
    -v GN="$L" '
    BEGIN{IGNORECASE=1}
    $0 ~ /^#/ {print; next}

    tolower($3)=="mrna" {
      sub(/ID=[^;]*/,"ID="NM,$9)
      if($9 ~ /Name=/){sub(/Name=[^;]*/,"Name="GN,$9)} else {$9=$9";Name="GN}
      print; next
    }

    tolower($3)=="cds" {
      sub(/Parent=[^;]*/,"Parent="NM,$9)
      sub(/ID=[^;]*/,"ID="NC,$9)
      if($9 ~ /Name=/){sub(/Name=[^;]*/,"Name="GN".CDS",$9)} else {$9=$9";Name="GN".CDS"}
      print; next
    }

    tolower($3)=="exon" {
      sub(/Parent=[^;]*/,"Parent="NM,$9)
      print; next
    }

    tolower($3)=="intron" {
      sub(/Parent=[^;]*/,"Parent="NM,$9)
      print; next
    }

    {print}
  ' "$OUT/exo.fixed.gff3" > "$OUT/${L}.gff3"

  echo "✅ Finished GFF3 for $L"
done

# ----------------------------
# 3) BUILD INTRON TABLE PER LOCUS
# ----------------------------
# This section reads CDS blocks from the GFF3 and infers intron intervals between them.
# It also carries phase information and calculates CodonIndexBefore (0-based-ish index of codons
# in the CDS before the intron boundary, based on summed CDS length).

for L in $LOCI; do
  OUT="$OUT_BASE/$L"
  GFF="$OUT/${L}.gff3"

  if [ ! -s "$GFF" ]; then
    echo "ERROR: missing cleaned GFF3 for $L ($GFF)"
    continue
  fi

  # Extract CDS rows and keep: Parent, start, end, strand, phase
  awk -F'\t' '
    tolower($3)=="cds"{
      if(match($9,/Parent=([^;]+)/,m)) print m[1]"\t"$4"\t"$5"\t"$7"\t"$8
    }' "$GFF" \
    | sort -k1,1 -k2,2n \
    > "$OUT/cds.sorted.tsv"

  # Determine strand from the first CDS row
  STRAND=$(awk 'NR==1{print $4}' "$OUT/cds.sorted.tsv")

  # Make a CDS list ordered in coding direction
  if [ "$STRAND" = "-" ]; then
    sort -k2,2nr "$OUT/cds.sorted.tsv" > "$OUT/cds.coding.tsv"
  else
    sort -k2,2n "$OUT/cds.sorted.tsv" > "$OUT/cds.coding.tsv"
  fi

  # (A) Basic intron coordinates using genomic order (mostly for debugging)
  awk '
    BEGIN{OFS="\t"}
    {
      p=$1; s=$2+0; e=$3+0; st=$4; np=$5+0
      if(p==pp){
        is=pe+1; ie=s-1
        if(ie>=is){
          i++
          il=ie-is+1
          iph=(3-np)%3
          print p,"intron_"i,is,ie,il,st,iph,np
        }
      } else { i=0 }
      pp=p; pe=e
    }' "$OUT/cds.sorted.tsv" > "$OUT/introns.basic.tsv"

  # (B) Intron table in coding order + CodonIndexBefore
  awk '
    BEGIN{OFS="\t"}
    {
      p=$1; s=$2+0; e=$3+0; st=$4; np=$5+0
      len=e-s+1

      if(NR==1){
        strand=st
      } else {
        # intron bounds depend on strand when walking in coding order
        if(strand=="+"){
          is=prev_e+1
          ie=s-1
        } else {
          is=e+1
          ie=prev_s-1
        }

        if(ie>=is){
          i++
          il=ie-is+1
          cb=int(prev_bp/3)       # codons completed before this intron
          iph=(3-np)%3            # intron phase relative to next CDS phase
          print p,"intron_"i,is,ie,il,strand,iph,np,cb
        }
      }

      prev_s=s
      prev_e=e
      prev_bp+=len
    }' "$OUT/cds.coding.tsv" > "$OUT/introns.with_codons.tsv"

  echo "✅ Built intron table for $L"
done

# ----------------------------
# 4) COMBINE ALL INTRONS INTO ONE TSV
# ----------------------------

OUTALL="$OUT_BASE/combined_introns.tsv"
mkdir -p "$(dirname "$OUTALL")"

# Header
echo -e "Gene\tParent\tIntron\tStart\tEnd\tLength\tStrand\tIntronPhase\tNextCDSPhase\tCodonIndexBefore" > "$OUTALL"

# Append per locus
for L in $LOCI; do
  f="$OUT_BASE/$L/introns.with_codons.tsv"
  if [ -s "$f" ]; then
    awk -v G="$L" 'BEGIN{OFS="\t"}{print G,$0}' "$f" >> "$OUTALL"
  fi
done

echo "✅ Combined introns written to: $OUTALL"

# ----------------------------
# 5) BUILD PEPTIDE MULTI-FASTA + RUN MAFFT
# ----------------------------
# This creates one peptide FASTA per locus with header "><LOCUS>",
# concatenates them, and runs MAFFT.

PEP_DIR="$MSA_OUT_DIR/peptides_by_locus"
IN_CAT="$MSA_OUT_DIR/peptides.all.fa"
OUT_ALN="$MSA_OUT_DIR/peptides.aln.fasta"

mkdir -p "$PEP_DIR"
rm -f "$PEP_DIR"/*.fa "$IN_CAT" "$OUT_ALN"

for L in $LOCI; do
  PEP="$IN_DIR/${L}_peptide_renamed.fasta"
  if [ ! -s "$PEP" ]; then
    echo "ERROR: missing peptide FASTA for $L ($PEP)"
    continue
  fi

  # Force a clean header that matches the locus ID
  awk -v H=">$L" 'NR==1{print H; next}{print}' "$PEP" > "$PEP_DIR/$L.fa"
done

if ls "$PEP_DIR"/*.fa >/dev/null 2>&1; then
  cat "$PEP_DIR"/*.fa > "$IN_CAT"
  mafft "${MAFFT_MODE[@]}" "$IN_CAT" > "$OUT_ALN"
  echo "✅ Aligned peptide FASTA written to: $OUT_ALN"
else
  echo "No peptide FASTAs found to align."
fi