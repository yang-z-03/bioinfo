
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

    score_interactions = True,     # optional: whether to score interactions or not. 
    iterations = 1000,             # denotes the number of shufflings performed in the analysis.
    threshold = 0.1,               # defines the min % of cells expressing a gene for this to be employed in the analysis.
    threads = 40,                  # number of threads to use in the analysis.
    debug_seed = 42,               # debug randome seed. To disable >=0.
    result_precision = 3,          # Sets the rounding for the mean values in significan_means.
    pvalue = 0.05,                 # P-value threshold to employ for significance.
    subsampling = False,           # To enable subsampling the data (geometri sketching).
    subsampling_log = False,       # (mandatory) enable subsampling log1p for non log-transformed data inputs.
    subsampling_num_pc = 100,      # Number of componets to subsample via geometric skectching (dafault: 100).
    subsampling_num_cells = 1000,  # Number of cells to subsample (integer) (default: 1/3 of the dataset).
    separator = '|',               # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    debug = False,                 # Saves all intermediate tables employed during the analysis in pkl format.
    
    output_path = '/home/yang-z/bioinfo/geo/GSE98638/cpdb',

    # Replaces the timestamp in the output files by a user defined string in the  (default: None).
    output_suffix = None
)