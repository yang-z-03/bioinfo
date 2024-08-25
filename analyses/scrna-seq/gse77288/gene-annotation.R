
library(tidyverse)
library(magrittr)
library(dplyr)
load("d:/projects/r/workspace.rdata")

library(AnnotationDbi)
library(GenomicFeatures)
library(EnsDb.Hsapiens.v86)

install_packages("EnsDb.Hsapiens.v86")

gene_table <- genes(EnsDb.Hsapiens.v86)

# > slotNames(gene_table)
# [1] "seqnames"        "ranges"          "strand"          "seqinfo"
# [5] "elementMetadata" "elementType"     "metadata"

# > ranges(gene_table)
# IRanges object with 63970 ranges and 0 metadata columns:
#                       start       end     width
#                   <integer> <integer> <integer>
#   ENSG00000223972     11869     14409      2541
#   ENSG00000227232     14404     29570     15167
#   ENSG00000278267     17369     17436        68
#   ENSG00000243485     29554     31109      1556
#   ENSG00000237613     34554     36081      1528
#               ...       ...       ...       ...
#   ENSG00000224240  26549425  26549743       319
#   ENSG00000227629  26586642  26591601      4960
#   ENSG00000237917  26594851  26634652     39802
#   ENSG00000231514  26626520  26627159       640
#   ENSG00000235857  56855244  56855488       245

# > strand(gene_table)
# factor-Rle of length 63970 with 31957 runs
#   Lengths: 1 2 1 1 3 2 1 5 2 6 1 1 1 2 1 3 ... 1 1 1 2 1 2 2 2 1 1 2 1 1 6 3 1
#   Values : + - + - + - + - + - + - + - + - ... - + - + - + - + - + - + - + - +
# Levels(3): + - *

# > seqinfo(gene_table)
# Seqinfo object with 357 sequences (1 circular) from GRCh38 genome:
#   seqnames seqlengths isCircular genome
#   1         248956422      FALSE GRCh38
#   10        133797422      FALSE GRCh38
#   11        135086622      FALSE GRCh38
#   12        133275309      FALSE GRCh38
#   13        114364328      FALSE GRCh38
#   ...             ...        ...    ...
#   LRG_741      231167      FALSE GRCh38
#   LRG_93        22459      FALSE GRCh38
#   MT            16569       TRUE GRCh38
#   X         156040895      FALSE GRCh38
#   Y          57227415      FALSE GRCh38

# > elementMetadata(gene_table)
# DataFrame with 63970 rows and 6 columns
#                         gene_id   gene_name           gene_biotype
#                     <character> <character>            <character>
# ENSG00000223972 ENSG00000223972     DDX11L1 transcribed_unproces..
# ENSG00000227232 ENSG00000227232      WASH7P unprocessed_pseudogene
# ENSG00000278267 ENSG00000278267   MIR6859-1                  miRNA
# ENSG00000243485 ENSG00000243485   MIR1302-2                lincRNA
# ENSG00000237613 ENSG00000237613     FAM138A                lincRNA
# ...                         ...         ...                    ...
# ENSG00000224240 ENSG00000224240     CYCSP49   processed_pseudogene
# ENSG00000227629 ENSG00000227629  SLC25A15P1 unprocessed_pseudogene
# ENSG00000237917 ENSG00000237917     PARP4P1 unprocessed_pseudogene
# ENSG00000231514 ENSG00000231514     FAM58CP   processed_pseudogene
# ENSG00000235857 ENSG00000235857     CTBP2P1   processed_pseudogene
#
#                 seq_coord_system      symbol                       entrezid
#                      <character> <character>                         <list>
# ENSG00000223972       chromosome     DDX11L1 100287596,100287102,727856,...
# ENSG00000227232       chromosome      WASH7P                             NA
# ENSG00000278267       chromosome   MIR6859-1                      102466751
# ENSG00000243485       chromosome   MIR1302-2            105376912,100302278
# ENSG00000237613       chromosome     FAM138A           654835,645520,641702
# ...                          ...         ...                            ...
# ENSG00000224240       chromosome     CYCSP49                             NA
# ENSG00000227629       chromosome  SLC25A15P1                             NA
# ENSG00000237917       chromosome     PARP4P1                             NA
# ENSG00000231514       chromosome     FAM58CP                             NA
# ENSG00000235857       chromosome     CTBP2P1                             NA
