# Intron mapping project

Comparative analysis of intron position conservation across CESA, CSLD, and related cellulose synthase–like genes in land plants and green algae.

This repository contains scripts, minimal input data, and documentation for extracting intron coordinates, mapping introns onto protein alignments, assessing positional conservation, and visualizing intron distributions across gene families.

Included workflows:

Exonerate-based intron extraction
A shell pipeline aligns peptide sequences to genomic loci, generates standardized GFF3 gene models, extracts intron coordinates, and builds peptide multiple-sequence alignments.

Intron-to-alignment mapping and matrix construction
An R script maps intron positions onto a protein multiple-sequence alignment and produces both tidy intron position tables and binary intron presence/absence matrices.

Intron presence heatmap visualization
An R script visualizes intron presence/absence patterns across genes as a heatmap from a binary intron matrix.

Permutation-based conservation testing
An R script evaluates intron positional conservation using permutation tests based on intron presence/absence matrices.

Alignment-based intron confidence scoring
An R script computes gap fraction and Shannon entropy around intron positions in the alignment and assigns confidence categories to intron mappings.