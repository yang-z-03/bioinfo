
from cellphonedb.utils import db_utils
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

# version of the databse
cpdb_version = 'v5'

# path where the input files to generate the database are located
cpdb_target_dir = '/home/yang-z/bioinfo/cpdb/cellphonedb.zip'

# check the input files
import anndata

adata = anndata.read_h5ad(
    '/home/yang-z/bioinfo/geo/GSE98638/cpdb/norm-counts.h5ad'
)

shape = adata.shape

cpdb_results = cpdb_statistical_analysis_method.call(

    # mandatory: CellphoneDB database zip file.
    cpdb_file_path = cpdb_target_dir,

    # mandatory: tsv file defining barcodes to cell label.
    meta_file_path = '/home/yang-z/bioinfo/geo/GSE98638/cpdb/metadata.tsv',

    # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
    counts_file_path = '/home/yang-z/bioinfo/geo/GSE98638/cpdb/norm-counts.h5ad',
    
    # defines the gene annotation in counts matrix.
    counts_data = 'hgnc_symbol',

    # optional: defines cell types and their active TFs.
    # active_tfs_file_path = active_tf_path,

    # optional (default: None): defines cells per microenvironment.
    # microenvs_file_path = microenvs_file_path,

    # optional: whether to score interactions or not. 
    score_interactions = True,

    # denotes the number of shufflings performed in the analysis.
    iterations = 1000,

    # defines the min % of cells expressing a gene for this to be employed in the analysis.
    # genes must exceed such expression to be considered.
    threshold = 0.1,

    threads = 40,
    debug_seed = 42,

    # Sets the rounding for the mean values in significan_means.
    result_precision = 3,
    
    # P-value threshold to employ for significance.
    pvalue = 0.05,

    # To enable subsampling the data (geometri sketching).
    subsampling = False,

    # (mandatory) enable subsampling log1p for non log-transformed data inputs.
    subsampling_log = False,

    # Number of componets to subsample via geometric skectching (dafault: 100).
    subsampling_num_pc = 100,

    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
    subsampling_num_cells = 1000,

    # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    separator = '|',
    
    # Saves all intermediate tables employed during the analysis in pkl format.
    debug = False,
    
    output_path = '/home/yang-z/bioinfo/geo/GSE98638/cpdb',

    # Replaces the timestamp in the output files by a user defined string in the  (default: None).
    output_suffix = None
)