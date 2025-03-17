
import os
import subprocess
from time import sleep
from printtbl import printtbl
from ansi import common_length
from ansi import info, warning, error
import scanpy as sc
import anndata as ad
import pickle

def hsize(value):
    units = ["B", "KiB", "MiB", "GiB", "TiB", "PiB"]
    size = 1024.0
    for i in range(len(units)):
        if (value / size) < 1:
            return "%.2f%s" % (value, units[i])
        value = value / size
    return 'too large'

import h5py
from scipy import sparse
from anndata.abc import CSRDataset as csr_dataset
from anndata.abc import CSCDataset as csc_dataset
from anndata.io import sparse_dataset
from anndata.io import read_elem, write_elem
from abc import ABC

class spsdataset(ABC):
    
    def __new__(cls, group):
        return sparse_dataset(group)
    
    @classmethod
    def __subclasshook__(cls, _c):
        return issubclass(_c, (csr_dataset, csc_dataset))

def read_metas(pth) -> ad.AnnData:
    
    attrs = ["obs", "var"]
    with h5py.File(pth) as f:
        adata = ad.AnnData(**{k: read_elem(f[k]) for k in attrs})
    return adata

def concat_on_disk(input_pths, output_pth):
    
    annotations = ad.concat([read_metas(pth) for pth in input_pths], join = 'outer', axis = 'obs')
    annotations.write_h5ad(output_pth)
    n_variables = annotations.shape[1]
    var_keys = annotations.var_names.tolist()
    annot_keys = annotations.obs_keys.tolist()
    obs_keys = []
    del annotations

    with h5py.File(output_pth, "a") as target:

        # read the data previously stored inside the raw.h5ad matrix.
        # and gives a adding point for later data to be added.
        
        tmpx = sparse.csr_matrix((0, n_variables), dtype = 'int16')
        tmpx.indptr = tmpx.indptr.astype('int64')
        write_elem(target, 'X', tmpx)
        mtx = spsdataset(target['X'])
        
        read_id = 0
        for p in input_pths:
            read_id += 1
            with h5py.File(p, 'r') as src:

                # ensure that empty columns (non-expressing genes) is inserted
                # accordingly to be appended to the matrix.
                # and well in the same order.
                
                blank_ad = ad.AnnData(shape = (0, 56834))
                blank_ad.var_names = var_keys
                blank_ad.X = sparse.csr_matrix((0, n_variables), dtype = 'int16')
                counts = ad.AnnData(read_elem(src['X']), var = read_elem(src['var']), obs = read_elem(src['obs']))
                reorder = ad.concat([blank_ad, counts], join = 'outer', axis = 'obs')
                ri16 = reorder.X.astype('int16')
                obs_keys += reorder.obs.index.tolist()
                
                print(
                    f'appending matrix for [{p}] ({read_id}/{len(input_pths)}), '
                    f'added {reorder.n_obs}, total {len(obs_keys)}'
                )
                
                assert np.sum(reorder.var.index == blank_ad.var.index) == \
                       len(reorder.var.index == blank_ad.var.index)
                
                mtx.append(ri16)
                del counts, reorder, ri16, blank_ad

    info('checking observation order ...')
    if len(obs_keys) != len(annot_keys):
        error('inconsistant sum of observations. routine will exit.')
    
    matched_obs = 0
    tested = 0
    for i in range(0, len(obs_keys), 1000):
        tested += 1
        if obs_keys[i] == annot_keys[i]: matched_obs += 1
    
    assert matched_obs == tested
    return


def display(registry, n_row, show_status = False):

    from qc import status_qc, styles_qc

    registry['qc'] = status_qc(registry)
    from list import styles_acctype, styles_number
    total_row = len(registry['qc'])
    nfrom = 0

    adatas = []
    n_total_obs = 0
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

        import pickle
        
        for s_acct, s_acc, s_name, s_qc in \
            zip(acctype, acc, nm, qcs):
            
            if s_qc == 'done':
                
                fpath = f'processed/qc/{s_name}.h5ad'
                szpath = f'processed/sizes/{s_name}.size'
                if os.path.exists(fpath):
                    
                    rprintc('reading ...')
                    n_obs = 0
                    n_vars = 0
                    if os.path.exists(szpath):
                        with open(szpath, 'rb') as szfile:
                            sz = pickle.load(szfile)
                            n_obs = sz['n_obs']
                            n_vars = sz['n_vars']
                    else:
                        a = sc.read_h5ad(f'processed/qc/{s_name}.h5ad')
                        n_obs = a.n_obs
                        n_vars = a.n_vars
                        with open(szpath, 'wb') as szfile:
                            pickle.dump({'n_obs': a.n_obs, 'n_vars': a.n_vars}, szfile)
                        
                        del a

                    should_append = False
                    if n_obs > 200 and n_vars > 200: 
                        
                        # a = sc.read_h5ad(f'processed/qc/{s_name}.h5ad')
                        # from scipy.sparse import csr_matrix
                        # if not isinstance(a.X, csr_matrix): a.X = csr_matrix(a.X)
                        # adatas.append(sc.read_h5ad(f'processed/qc/{s_name}.h5ad'))
                        
                        adatas.append(fpath)
                        should_append = True
                        n_total_obs += n_obs
                        
                    fprintc(f'[obs] \033[32;1m{n_obs}\033[0m * [var] \033[31;1m{n_vars}\033[0m ' + 
                            ('*' if should_append else ''), l = 62)
                
                else: fprintc('cannot find file.')
                
            else: fprintc('skipped.')

        nfrom += n_row
        print('\r')

    with open('processed/integrated/qc-datasets.pkl', 'wb') as qcd:
        pickle.dump(adatas, qcd)
    
    print(f'merged data has [obs] \033[32;1m{n_total_obs}\033[0m observations.')
    if not os.path.exists('processed/integrated/counts.h5ad'):

        info('begin merging adatas ...')

        # the internal experimental function for concat on disk stopped at somewhat
        # 150gb file size. and i am unable to debug that ...
        
        ad.experimental.concat_on_disk(
            in_files = adatas, max_loaded_elems = 5e9, # about 20 gb usage (5g elementes at the same time.)
            out_file = 'processed/integrated/counts.h5ad',
            axis = 'obs', # concatentate observations
            join = 'outer', # for sparse matrix, 0's will be filled for variables that may not exist
            merge = None, uns_merge = None, # throw var metadata for speed. the metadata of genes is easily added later.
            label = None, # there is already a 'accession' column, no need here.
            fill_value = None,
            pairwise = None
        )

        # almost sure this will run out of memory.
        
        # a = ad.concat(adatas, axis = 'obs', join = 'outer')
        # a.write_h5ad('processed/integrated/raw.h5ad')

        # let's try the hacky part ourself :(
        
        concat_on_disk(adatas, 'processed/integrated/counts.h5ad')
        hrsize = hsize(os.path.getsize('processed/integrated/counts.h5ad'))
        info(f'saved to processed/integrated/raw.h5ad. [{hrsize}]')

    else:
        
        hrsize = hsize(os.path.getsize('processed/integrated/counts.h5ad'))
        info(f'processed/integrated/raw.h5ad exists. [{hrsize}]')

    return # ends here.

    a = None
    
    import pandas as pd
    gname = a.var_names.tolist()
    gtable = pd.read_table('genome.tsv')
    finder = gtable[['gene']].values.transpose()[0].tolist()
    gind = ['g' + str(i + 1) for i in range(len(finder))]
    gtable['.ugene'] = gind
    gtable = gtable.set_index('.ugene')
    
    full.var['.seqid'] = gtable['.seqid']
    full.var['.start'] = gtable['.start']
    full.var['.end'] = gtable['.end']
    full.var['.strand'] = gtable['.strand']
    full.var['.name'] = gtable['gene']
    full.var['ensembl'] = gtable['ensembl']
    full.var['refseq'] = gtable['refseq']
    full.var['uniprot'] = gtable['uniprot']
    full.var['mgi'] = gtable['.mgi']

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
    a.var['hvg.all'] = hvga['hvg']
    a.var['hvg.all'] = a.var['hvg.all'].fillna(False)
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

    pass
