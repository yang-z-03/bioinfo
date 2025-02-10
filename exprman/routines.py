
import os
from io import TextIOWrapper as textio
import omicverse as ov
from scipy.sparse import csr_matrix
import pandas as pd
import scanpy as sc
import numpy as np

gtable = pd.read_table('genome.tsv')

def compile_mtx(src, prefix, props, clog: textio):

    # scanpy can only recognize mtx in non-gzipped format. so we should
    # decompress them if possible.
    
    from sh import gunzip

    if os.path.exists(f'{src}/{prefix}matrix.mtx.gz'):
        gunzip(f'{src}/{prefix}matrix.mtx.gz')
    if os.path.exists(f'{src}/{prefix}genes.tsv.gz'):
        gunzip(f'{src}/{prefix}genes.tsv.gz')
    if os.path.exists(f'{src}/{prefix}barcodes.tsv.gz'):
        gunzip(f'{src}/{prefix}barcodes.tsv.gz')

    adata = sc.read_10x_mtx(
        src, var_names = 'gene_ids',
        gex_only = True, make_unique = True,
        prefix = prefix
    )

    final = process_matrix(props, clog, adata)
    
    del adata
    return final


def process_matrix(props, clog, adata, force_filter = False):

    clog.writelines([f'    [i] raw matrix read: {adata.n_obs} obs * {adata.n_vars} var \n'])

    # if more than 50000 cells in a single matrix, we just believe that it is
    # an unfiltered raw matrix. we should roughly filter the empty droplets
    if adata.n_obs > 50000 or force_filter:
        valid = (adata.X.sum(axis = 1).transpose() > 200).tolist()[0]
        adata_f = adata[valid, :]
    else: adata_f = adata

    clog.writelines([f'    [i] filtered matrix: {adata_f.n_obs} obs * {adata_f.n_vars} var \n'])

    # append the sample name to the barcode
    adata_f.obs_names = props['name'] + ':' + adata_f.obs_names

    # map gene naming
    gname = adata_f.var_names.tolist()
    names = []
    gmask = []

    # here, we just add another condition to test whether the gname list is appropriate
    # ensembl format. if it is not, we try to map genes directly onto the names.
    # though i specify the var_names should be 'gene_ids', it may occur exceptions
    # where there are man-made references containing two or more species. by convention
    # in these double species reference, the 'gene_ids' should be 'mm10_ENSMUSG...'
    # or just name of the genes.

    do_start_with_mm10 = 0
    do_ensembl = 0
    for gn in gname:
        if gn.startswith('ENSMUSG'): do_ensembl += 1
    for gn in gname:
        if '_ENSMUSG' in gn: do_start_with_mm10 += 1
    
    if do_start_with_mm10 > 0:
        gname_copy = []
        for gn in gname:
            if '_' in gn: gname_copy += [gn.split('_')[1]]
            else: gname_copy += [gn]
        gname = gname_copy
    
    finder = gtable[['ensembl' if (do_start_with_mm10 + do_ensembl > 0) else 'gene']] \
        .values.transpose()[0].tolist()
    not_in_list = []

    alias_table = {}
    # make the name alias table.
    if (do_start_with_mm10 + do_ensembl) == 0:

        finder_alias = gtable[['alias']].values.transpose()[0].tolist()
        duplicates = []

        for x in range(len(finder_alias)):
            alias = finder_alias[x]
            if not isinstance(alias, str): continue
            if len(alias) > 0:
                spl = alias.split(';')
                for y in spl:
                    if y not in alias_table.keys():
                        alias_table[y] = x
                    elif y not in duplicates:
                        duplicates += [y]
        
        for x in duplicates:
            del alias_table[x]

    for x in gname:
        if x in finder:
            gmask += [True]
            names += ['g' + str(finder.index(x) + 1)]
        elif x in alias_table.keys(): # if using ensembl, this is an empty table.
            gmask += [True]
            names += ['g' + str(alias_table[x] + 1)]
        else:
            gmask += [False]
            not_in_list += [x]
    

    clog.writelines([
        f'    [i] {len(not_in_list)} genes (out of {len(gname)}) not in the reference gene list ({len(finder)}) \n',
        f'    [i] total {len(names)} genes mapped. {len(np.unique(names))} unique genes. \n'
    ])

    final = adata_f[:, gmask].copy()
    final.var_names = names
    # remove duplicated genes
    final = final[:, final.var_names.duplicated() == False]

    # attach cell metadata onto the obs slot.
    final.obs['sample'] = props['name']
    final.obs['accession'] = props['accession']
    final.obs['dataset'] = props['dataset']
    final.obs['tissue'] = props['tissue']
    final.obs['genetic'] = props['genetic']
    final.obs['strain'] = props['strain']
    final.obs['sort'] = props['sort']
    final.obs['expm'] = props['expm']
    final.obs['age'] = props['age']
    final.obs['sex'] = props['sex']
    final.obs['ms'] = props['ms']
    final.obs['fs'] = props['fs']
    final.obs['barcode'] = final.obs_names
    final.obs['ubc'] = [props['accession'] + ':' + str(x + 1) for x in range(final.n_obs)]
    final.obs_names = [props['accession'] + ':' + str(x + 1) for x in range(final.n_obs)]

    if not isinstance(final.X, csr_matrix):
        final.X = csr_matrix(final.X)

    del adata_f
    return final


def compile_mtxz(src, prefix, props, clog: textio):
    
    adata = sc.read_10x_mtx(
        src, var_names = 'gene_ids',
        gex_only = True, make_unique = True,
        prefix = prefix
    )

    final = process_matrix(props, clog, adata)
    
    del adata
    return final

def compile_h5r(src, h5name, props, clog: textio):

    adata = sc.read_10x_h5(
        f'{src}/{h5name}',
        gex_only = True
    )

    lst = adata.var['gene_ids'].tolist()
    lsstab = [y.split('.')[0] for y in lst]
    adata.var_names = lsstab
    del adata.var

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    final = process_matrix(props, clog, adata, force_filter = True)
    
    del adata
    return final


def compile_h5f(src, h5name, props, clog: textio):
    
    adata = sc.read_10x_h5(
        f'{src}/{h5name}',
        gex_only = True
    )

    lst = adata.var['gene_ids'].tolist()
    lsstab = [y.split('.')[0] for y in lst]
    adata.var_names = lsstab
    del adata.var
    
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    final = process_matrix(props, clog, adata, force_filter = False)
    
    del adata
    return final