
import os
import subprocess
from time import sleep
from printtbl import printtbl
from ansi import common_length
from ansi import info, warning, error
import omicverse as ov
import scanpy as sc
import anndata as ad

def hsize(value):
    units = ["B", "KiB", "MiB", "GiB", "TiB", "PiB"]
    size = 1024.0
    for i in range(len(units)):
        if (value / size) < 1:
            return "%.2f%s" % (value, units[i])
        value = value / size
    return 'too large'

def display(registry, n_row, show_status = False):

    from qc import status_qc, styles_qc

    registry['qc'] = status_qc(registry)
    from list import styles_acctype, styles_number
    total_row = len(registry['qc'])
    nfrom = 0

    adatas = []

    while nfrom < total_row:

        cols = printtbl(
            registry, styles = { 
                'type': styles_acctype, 'accession': styles_number,
                'qc': styles_qc
            },
            appends = { 
                'type': 9, 'accession': 11,
                'qc': 11
            }, nfrom = nfrom, nto = nfrom + n_row - 1
        )

        # move the cursor to the first item.
        from ansi import ansi_move_cursor
        ansi_move_cursor(-(min(total_row, nfrom + n_row) - nfrom), cols)

        def print_final(x, reallen):
            print(x, end = '')
            ansi_move_cursor(1, -reallen)

        def print_replace(x, reallen):
            print(x, end = '')
            ansi_move_cursor(0, -reallen)
        
        def fprintc(x, l = 40, reallen = 40): print_final(common_length(x, l), reallen)
        def rprintc(x, l = 40, reallen = 40): print_replace(common_length(x, l), reallen)

        qcs = registry['qc'][nfrom : min(total_row, nfrom + n_row)]
        acctype = registry['type'][nfrom : min(total_row, nfrom + n_row)]
        acc = registry['accession'][nfrom : min(total_row, nfrom + n_row)]
        nm = registry['name'][nfrom : min(total_row, nfrom + n_row)]
        
        for s_acct, s_acc, s_name, s_qc in \
            zip(acctype, acc, nm, qcs):
            
            if s_qc == 'done':
                rprintc('reading ...')
                a = sc.read_h5ad(f'processed/qc/{s_name}.h5ad')
                adatas += [a]
                fprintc(f'[obs] \033[32;1m{a.n_obs}\033[0m * [var] \033[31;1m{a.n_vars}\033[0m', l = 62)
            
            else: fprintc('skipped.')

        nfrom += n_row
        print('\r')
    
    info('begin merging adatas ...')
    a = ad.concat(adatas, join = 'outer')
    print(f'merged adata: [obs] \033[32;1m{a.n_obs}\033[0m * [var] \033[31;1m{a.n_vars}\033[0m')

    import pandas as pd
    gname = a.var_names.tolist()
    gtable = pd.read_table('genome.tsv')
    finder = gtable[['gene']].values.transpose()[0].tolist()
    gind = ['g' + str(i + 1) for i in range(len(finder))]
    gtable['.ugene'] = gind
    gtable = gtable.set_index('.ugene')

    a.var['.seqid'] = gtable['.seqid']
    a.var['.start'] = gtable['.start']
    a.var['.end'] = gtable['.end']
    a.var['.strand'] = gtable['.strand']
    a.var['.name'] = gtable['gene']
    a.var['ensembl'] = gtable['ensembl']
    a.var['refseq'] = gtable['refseq']
    a.var['uniprot'] = gtable['uniprot']
    a.var['mgi'] = gtable['.mgi']

    # here, we rewrite the batch key.
    a.obs['batch'] = a.obs['accession'].str.cat(a.obs['sample'], sep = ':').astype('category')

    print('')

    warning('reading all hvgs ...')
    # reading all hvgs from all samples:
    import pickle
    hvgfs = os.listdir('processed/hvg')
    hvglists = []
    for fn in hvgfs:
        with open('processed/hvg/' + fn, 'rb') as f:
            hvglists += [pickle.load(f)]

    hvgall = []
    for x in hvglists: hvgall += list(x)
    hvgall = list(set(hvgall))
    info(f'you have loaded {len(hvglists)} hvg sets, a total of {len(hvgall)} variable genes.')
    
    # attaching hvg.all
    hvga = { '.ugene': hvgall, 'hvg': [True for i in range(len(hvgall))]}
    hvga = pd.DataFrame(hvga)
    hvga = hvga.set_index('.ugene')
    a.var['hvg.all'] = False
    a.var['hvg.all'] = hvga['hvg']
    a.var['hvg.all'] = a.var['hvg.all'].astype('category')
    info(f'all non-redundant variable genes ({len(hvgall)}) is attached to var[hvg.all]')

    print('')

    # rerun top 4000 genes for dimension reduction.

    # singularity value errors may happen when running with the default parameter
    # span = 0.3 using seurat_v3 hvg algorithm when fitting the loess regression.

    # the problem can occur when batch key has many cells in each batch. 
    # Increasing the span from the default of 0.3 to 0.5 seems to have "fixed" 
    # the error. Increasing the filtering stringency for lowly expressed genes 
    # (to min_gene = 500, min_cells = 10) also gets rid of the error.

    warning('rerunning hvg selection for top 4000 hvgs for accession ...')
    a.layers['counts'] = a.X
    seurat_v3_hvg = sc.pp.highly_variable_genes(
        a, flavor = 'seurat_v3', layer = 'counts',
        n_top_genes = 4000, subset = False, inplace = False,
        batch_key = 'accession', span = 0.5
    )

    del a.layers
    a.var['hvg.4000'] = seurat_v3_hvg['highly_variable']
    a.uns['hvg.4000'] = seurat_v3_hvg
    info('done.')

    print('')
    warning('atlas variable table: \n')
    print(a.var[['.seqid', '.start', '.end', '.strand', '.name', 'ensembl']])

    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_colwidth', 90)
    print('')
    warning('atlas observation table: \n')
    print(a.obs)

    from scipy.sparse import csr_matrix
    if not isinstance(a.X, csr_matrix): a.X = csr_matrix(a.X)
    
    print('')
    warning('writing to disk (raw.h5ad) ...')
    a.write_h5ad('processed/integrated/raw.h5ad')
    hrsize = hsize(os.path.getsize('processed/integrated/raw.h5ad'))
    info(f'saved to processed/integrated/raw.h5ad. [{hrsize}]')

    # print('')
    # warning('writing raw count matrix only (counts.h5ad) ...')
    # counts = ad.AnnData(a.X.copy())
    # counts.obs_names = a.obs_names
    # counts.var_names = a.var_names
    # counts.obs['batch'] = a.obs['batch'].astype('category')
    # counts.write_h5ad('processed/integrated/counts.h5ad')
    # hrsize = hsize(os.path.getsize('processed/integrated/counts.h5ad'))
    # info(f'saved to processed/integrated/counts.h5ad. [{hrsize}]')
    
    pass
