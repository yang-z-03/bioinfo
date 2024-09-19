
# cellphonedb is a suite to study cell-cell communication from single-cell
# transcriptomics data. identifying ligand–receptor interactions from scrna-seq
# requires both the annotation of complex ligand–receptor relationships from
# the literature (i.e. database) and a statistical method that integrates the
# resource with scrna-seq data and selects relevant interactions from the
# dataset (ie.e tool). cellphonedb is composed of two units, a database
# and a tool.

# cellphonedb database is a publicly available repository of curated receptors,
# ligands and their interactions. subunit architecture is included for both
# ligands and receptors, representing heteromeric complexes accurately. this
# is crucial, as cell-cell communication relies on multi-subunit protein
# complexes that go beyond the binary representation used in most databases and
# studies. it integrates new manually reviewed interactions with evidenced
# roles in cell-cell communication together with existing datasets that pertain
# to cellular communication (such as shilts et al. 2022 and
# kanemaru et al. 2023). recently, the database expanded to include non-protein
# molecules acting as ligands.

# the database can be used to search for a particular ligand/receptor or in
# combination with the tool to interrogate your own single-cell
# transcriptomics data.


# 1.  INPUTS ===================================================================

# The statistial method accepts 4 input files (3 mandatory).
#
#   cpdb_file_path: (mandatory) path to the database cellphonedb.zip.
#   meta_file_path: (mandatory) path to the meta file linking cell barcodes
#       to cluster labels metadata.tsv.
#   counts_file_path: (mandatory) paths to normalized counts file
#       (not z-transformed), either in text format or h5ad. (log-norm counts)
#   microenvs_file_path (optional) path to microenvironment file that groups
#       cell types/clusters by microenvironments. When providing a
#       microenvironment file, CellphoneDB will restrict the interactions to
#       those cells within the microenvironment.
