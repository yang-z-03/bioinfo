
# gene2ensembl                                    recalculated daily
# -----------------------------------------------------------------------------
#            This file reports matches between NCBI and Ensembl annotation
#            based on comparison of rna and protein features.
#
#            Matches are collected as follows.
#            For a protein to be identified as a match between RefSeq and
#            Ensembl, there must be at least 80% overlap between the two.
#            Furthermore, splice site matches must meet certain conditions:
#            either 60% or more of the splice sites must match, or there may
#            be at most one splice site mismatch.
#
#            For rna features, the best match between RefSeq and Ensembl is
#            selected based on splice site and overlap comparisons. For coding
#            transcripts, there is no minimum threshold for reporting other than
#            the protein comparison criteria above. For non-coding transcripts,
#            the splice site criteria are the same as for protein matching, but
#            the overlap threshold is reduced to 50%.
#
#            Furthermore, both the rna and the protein features must meet these
#            minimum matching criteria to be considered a good match.  In
#            addition, only the best matches will be reported in this file.
#            Other matches that satisified the matching criteria but were
#            not the best matches will not be reported in this file.
#
#            A summary report of species that have been compared is contained
#            in another FTP file, README_ensembl (see below).
#
#            Ensembl gene identifiers are also reported in the dbXrefs column
#            in the gene_info FTP file.  Due to differences in how these files
#            are processed, the Ensembl gene identifiers in these two files
#            may not be in complete concordance.
#
#            More notes about this file:
#
#            tab-delimited
#            one line per match between RefSeq and Ensembl rna/protein
#            Column header line is the first line in the file.
#
#
# -----------------------------------------------------------------------------
#
# tax_id:
#            the unique identifier provided by NCBI Taxonomy
#            for the species or strain/isolate
#
# GeneID:
#            the unique identifier for a gene
#
# Ensembl_gene_identifier:
#            the matching Ensembl identifier for the gene
#
# RNA nucleotide accession.version:
#            the identifier for the matching RefSeq rna
#            will be null (-) if only the protein matched
#
# Ensembl_rna_identifier:
#            the identifier for the matching Ensembl rna
#            may include a version number
#            will be null (-) if only the protein matched
#
# protein accession.version:
#            the identifier for the matching RefSeq protein
#            will be null (-) if only the mRNA matched
#
# Ensembl_protein_identifier:
#            the identifier for the matching Ensembl protein
#            may include a version number
#            will be null (-) if only the mRNA matched

# Homo sapiens Annotation Release 110	GRCh38.p14	112	GRCh38.p14

ncbi_ensembl_accession <-
  function(accession_file = "refseq/convertions/ensembl.gz") {

    cnames <- c("tax", "refseq", "gene", "version_rna", "rna",
                "version_protein", "protein")
    table_data <- utils::read.delim(accession_file, header = FALSE, sep = "\t",
                                    comment.char = "#", quote = "",
                                    blank.lines.skip = TRUE) |>
      `colnames<-`(cnames) |>
      as.data.table()

    return(table_data)
  }

ncbi_ensembl_accession_human <- function(g2e_table) {
  g2e_table <- g2e_table[tax == 9606]
}
