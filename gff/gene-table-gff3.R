
# this is the feature table for refseq 9606

# > pairs <- gff3_pairs('refseq/9606/9606_genomic_refseq.gff.gz')
# > pairs
#                                         element             parent       n
# 1                                           CDS     C_gene_segment     158
# 2                                           CDS     D_gene_segment      61
# 3                                           CDS     J_gene_segment     125
# 4                                           CDS     V_gene_segment    1224
# 5                                           CDS               mRNA 1836612
# 6                                C_gene_segment               gene      37
# 7                                C_gene_segment         pseudogene       7
# 8                                D_gene_segment               gene      61
# 9                                J_gene_segment               gene     117
# 10                               J_gene_segment         pseudogene      11
# 11                                RNase_MRP_RNA               gene       1
# 12                                  RNase_P_RNA               gene       2
# 13                               V_gene_segment               gene     365
# 14                               V_gene_segment         pseudogene     299
# 15                                        Y_RNA               gene       4
# 16                                antisense_RNA               gene      39
# 17                                     enhancer               gene       6
# 18                                         exon     C_gene_segment     165
# 19                                         exon     D_gene_segment      61
# 20                                         exon     J_gene_segment     128
# 21                                         exon      RNase_MRP_RNA       1
# 22                                         exon        RNase_P_RNA       2
# 23                                         exon     V_gene_segment    1224
# 24                                         exon              Y_RNA       4
# 25                                         exon      antisense_RNA     138
# 26                                         exon            lnc_RNA  120466
# 27                                         exon               mRNA 2003634
# 28                                         exon              miRNA    3218
# 29                                         exon              ncRNA      54
# 30                                         exon primary_transcript    2139
# 31                                         exon         pseudogene    9614
# 32                                         exon               rRNA      80
# 33                                         exon              scRNA       4
# 34                                         exon              snRNA     172
# 35                                         exon             snoRNA    1300
# 36                                         exon               tRNA     734
# 37                                         exon     telomerase_RNA       1
# 38                                         exon         transcript  164636
# 39                                         exon          vault_RNA       4
# 40                                      lnc_RNA               gene   32792
# 41                                      lnc_RNA         pseudogene      52
# 42                                         mRNA               gene  144481
# 43                                        miRNA primary_transcript    3218
# 44                                        ncRNA               gene      54
# 45                           primary_transcript               gene    2139
# 46                                         rRNA               gene      80
# 47                                        scRNA               gene       4
# 48                                        snRNA               gene     172
# 49                                       snoRNA               gene    1300
# 50                                         tRNA               gene     691
# 51                               telomerase_RNA               gene       1
# 52                                   transcript               gene   12833
# 53                                   transcript         pseudogene    2141
# 54                                    vault_RNA               gene       4
# 55                                  CAAT_signal                          6
# 56                                 CAGE_cluster                         96
# 57                   DNaseI_hypersensitive_site                        177
# 58                                       D_loop                          1
# 59                      GC_rich_promoter_region                         15
# 60                                     TATA_box                         30
# 61                            biological_region                     137091
# 62                                   cDNA_match                      26404
# 63                                   centromere                         24
# 64                        chromosome_breakpoint                         14
# 65                             conserved_region                        138
# 66                                direct_repeat                         20
# 67                             dispersed_repeat                          5
# 68                                     enhancer                     115157
# 69                    enhancer_blocking_element                         59
# 70               epigenetically_modified_region                         12
# 71                                         gene                      48104 *
# 72                    imprinting_control_region                          2
# 73                                    insulator                         12
# 74                         locus_control_region                         14
# 75                                        match                     161162
# 76                       matrix_attachment_site                         18
# 77                 meiotic_recombination_region                        368
# 78                               microsatellite                          5
# 79                                minisatellite                         13
# 80                 mitotic_recombination_region                         54
# 81                       mobile_genetic_element                        204
# 82  non_allelic_homologous_recombination_region                        514
# 83                     nucleotide_cleavage_site                          4
# 84                             nucleotide_motif                        623
# 85                        origin_of_replication                         87
# 86                                     promoter                        438
# 87                         protein_binding_site                       1418
# 88                                   pseudogene                      19252
# 89                        recombination_feature                        632
# 90                                       region                        709
# 91                            regulatory_region                          2
# 92                    repeat_instability_region                         62
# 93                                repeat_region                         10
# 94                replication_regulatory_region                         11
# 95                       replication_start_site                          3
# 96                             response_element                         20
# 97                          sequence_alteration                         30
# 98                 sequence_alteration_artifact                         10
# 99                          sequence_comparison                          2
# 100                            sequence_feature                       1964
# 101                sequence_secondary_structure                         14
# 102                                    silencer                      35932
# 103                               tandem_repeat                         64
# 104       transcriptional_cis_regulatory_region                       2259


gff3_list_regions_refseq <- function(gff_table) {
  gff3_list_attributes(gff_table, filter_feature = "region",
                       patterns = c("Dbxref", "Name", "chromosome", "gbkey",
                                    "genome", "mol_type"),
                       column_names = c("xref", "name", "chromosome",
                                        "genbank", "genome", "molecular_type"))
}

gff3_list_chromosomes_refseq <- function(gff_table) {
  gff3_list_regions_refseq(gff_table)[genome == "chromosome"] # nolint
}

gff3_list_genes_refseq <- function(gff_table) {
  gff3_list_attributes(gff_table, filter_feature = "gene",
                       patterns = c("ID", "Dbxref", "Name", "description",
                                    "gbkey", "gene", "gene_biotype"),
                       column_names = c("id", "xref", "name", "description",
                                        "genbank", "gene", "biotype"))
}

gff3_list_genes_ensembl <- function(gff_table) {
  gff3_list_attributes(gff_table, filter_feature = "gene",
                       patterns = c("ID", "Name", "biotype", "description",
                                    "gene_id", "logic_name", "version"),
                       column_names = c("id", "name", "biotype", "description",
                                        "gene_id", "logic_name", "version"))
}

gff3_list_chromosomes_ensembl <- function(gff_table) {
  gff3_list_attributes(gff_table, filter_feature = "chromosome",
                       patterns = c("ID", "Alias"),
                       column_names = c("id", "xref"))
}