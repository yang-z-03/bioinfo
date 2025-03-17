
import os
from io import TextIOWrapper as textio
import omicverse as ov
from scipy.sparse import csr_matrix
import pandas as pd
import scanpy as sc
import numpy as np
import pickle

gtable = pd.read_table('genome.tsv')
REFINE_FINDER = False

# we will construct the gene table for ensembl id and gene name alias at once during
# the startup process of this script. then we will be able to directly calls for
# its ugene id in all the next times.

if os.path.exists('configs/ensembl-finder.pkl') and os.path.exists('configs/name-finder.pkl'):
    with open('configs/ensembl-finder.pkl', 'rb') as fens:
        _ensembl_finder = pickle.load(fens)
    with open('configs/name-finder.pkl', 'rb') as fname:
        _name_finder = pickle.load(fname)
        
    print(f'Reading ensembl finder dictionary ... {len(_ensembl_finder)}')
    print(f'Reading gene symbol finder dictionary ... {len(_name_finder)}')
    
else:
    _ensembl_list = gtable[['ensembl']].values.transpose()[0].tolist()
    _ensembl_finder = {}
    
    _name_list = gtable[['gene']].values.transpose()[0].tolist()
    _name_finder = {}
    
    _alias_list = gtable[['alias']].values.transpose()[0].tolist()
    _duplicates = []
    _alias_finder = {}
    for x in range(len(_alias_list)):
        alias = _alias_list[x]
        if not isinstance(alias, str): continue
        if len(alias) > 0:
            spl = alias.split(';')
            for y in spl:
                if y not in _alias_finder.keys():
                    _alias_finder[y] = x
                elif y not in _duplicates:
                    _duplicates += [y]
    
    for x in _duplicates:
        del _alias_finder[x]
    
    for k in _alias_finder.keys():
        _name_finder[k] = 'g' + str(1 + _alias_finder[k])
    
    for i in range(len(_name_list)):
        _name_finder[_name_list[i]] = 'g' + str(1 + i)
    
    for i in range(len(_ensembl_list)):
        _ensembl_finder[_ensembl_list[i]] = 'g' + str(1 + i)
        
    print(f'Constructing ensembl finder dictionary ... {len(_ensembl_finder)}')
    print(f'Constructing gene symbol finder dictionary ... {len(_name_finder)}')

print()

def compile_mtx(src, prefix, props, clog: textio, raw = False):

    # scanpy can only recognize mtx in non-gzipped format. so we should
    # decompress them if possible.

    gpath = f'{src}/{prefix}genes.tsv.gz' if os.path.exists(f'{src}/{prefix}genes.tsv.gz') \
            else f'{src}/{prefix}genes.tsv'
    
    if not os.path.exists(gpath):
        return None, ['mismatch between assumed data type and actual file?', gpath]
    
    gfile = pd.read_table(gpath, header = None)
    
    if len(gfile.columns) == 3:
        clog.writelines([
            '    [i] there are 3 columns in the genes.tsv.gz\n',
            '        [1] ' + str(gfile[0].tolist()[0:5]) + '\n',
            '        [2] ' + str(gfile[1].tolist()[0:5]) + '\n',
            '        [3] ' + str(gfile[2].tolist()[0:5]) + '\n',
            '        redirecting from [mtx] to [mtxz] ... \n',
        ])
        os.rename(gpath, gpath.replace('genes.tsv', 'features.tsv'))
        return compile_mtxz(src, prefix, props, clog, raw)
    
    elif len(gfile.columns) == 2:
        
        # check the format, the first col should be ENSMUSG index, and the
        # second should be gene name.
        
        clog.writelines([
            '    [i] there are 2 columns in the genes.tsv.gz\n',
            '        [1] ' + str(gfile[0].tolist()[0:5]) + '\n',
            '        [2] ' + str(gfile[1].tolist()[0:5]) + '\n'
        ])

        # refine our finders:

        if REFINE_FINDER:
            ensembls = gfile[0].tolist()
            names = gfile[0].tolist()
            for ens, nm in zip(ensembls, names):
                if (nm not in _name_finder.keys()) and (ens in _ensembl_finder.keys()):
                    _name_finder[nm] = _ensembl_finder[ens]
    
            with open('configs/ensembl-finder.pkl', 'wb') as fens:
                pickle.dump(_ensembl_finder, fens)
            with open('configs/name-finder.pkl', 'wb') as fname:
                pickle.dump(_name_finder, fname)
            
        pass
    
    else:
        return None, [
            f'illegal genes.tsv', 
            f'{len(gfile.columns)} columns is observed.'
        ] + [f'[{_x + 1}] ' + str(gfile[_x].tolist()[0:5]) for _x in range(len(gfile.columns))]

    from sh import gunzip

    if os.path.exists(f'{src}/{prefix}matrix.mtx.gz'):
        gunzip(f'{src}/{prefix}matrix.mtx.gz')
    if os.path.exists(f'{src}/{prefix}genes.tsv.gz'):
        gunzip(f'{src}/{prefix}genes.tsv.gz')
    if os.path.exists(f'{src}/{prefix}barcodes.tsv.gz'):
        gunzip(f'{src}/{prefix}barcodes.tsv.gz')

    try:
        adata = sc.read_10x_mtx(
            src, var_names = 'gene_ids',
            gex_only = True, make_unique = True,
            prefix = prefix
        )
    except:
        return None, ['reading mtx error.']

    final = process_matrix(props, clog, adata, force_filter = raw)
    
    del adata
    return final, []


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

    do_start_with_mm10 = False
    do_ensembl = False
    for gn in gname:
        if gn.startswith('ENSMUSG'):
            do_ensembl = True
            continue
    
    for gn in gname:
        if '_ENSMUSG' in gn: 
            do_start_with_mm10 += True
            continue
    
    if do_start_with_mm10 > 0:
        gname_copy = []
        for gn in gname:
            if '_' in gn: gname_copy += [gn.split('_')[1]]
            else: gname_copy += [gn]
        gname = gname_copy
    
    not_in_list = []

    finder = {}
    if do_start_with_mm10 or do_ensembl:
        finder = _ensembl_finder
    else: finder = _name_finder

    keys = finder.keys()
    for x in gname:
        if x in keys:
            gmask.append(True)
            names.append(finder[x])
        else:
            gmask.append(False)
            not_in_list.append(x)
    
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


def compile_mtxz(src, prefix, props, clog: textio, raw = False):

    gpath = f'{src}/{prefix}features.tsv.gz' if os.path.exists(f'{src}/{prefix}features.tsv.gz') \
            else f'{src}/{prefix}features.tsv'
    if not os.path.exists(gpath):
        return None, ['mismatch between assumed data type and actual file?', gpath]
    gfile = pd.read_table(gpath, header = None)
    
    if len(gfile.columns) == 3:
        clog.writelines([
            '    [i] there are 3 columns in the features.tsv\n',
            '        [1] ' + str(gfile[0].tolist()[0:5]) + '\n',
            '        [2] ' + str(gfile[1].tolist()[0:5]) + '\n',
            '        [3] ' + str(gfile[2].tolist()[0:5]) + '\n',
        ])

        # refine our finders:
        if REFINE_FINDER:
            ensembls = gfile[0].tolist()
            names = gfile[0].tolist()
            for ens, nm in zip(ensembls, names):
                if (nm not in _name_finder.keys()) and (ens in _ensembl_finder.keys()):
                    _name_finder[nm] = _ensembl_finder[ens]
            
            with open('configs/ensembl-finder.pkl', 'wb') as fens:
                pickle.dump(_ensembl_finder, fens)
            with open('configs/name-finder.pkl', 'wb') as fname:
                pickle.dump(_name_finder, fname)
        
        pass

    elif len(gfile.columns) == 2:
        
        # check the format, the first col should be ENSMUSG index, and the
        # second should be gene name.
        
        clog.writelines([
            '    [i] there are 2 columns in the features.tsv\n',
            '        [1] ' + str(gfile[0].tolist()[0:5]) + '\n',
            '        [2] ' + str(gfile[1].tolist()[0:5]) + '\n',
            '        redirecting from [mtxz] to [mtx] ... \n',
        ])
        os.rename(gpath, gpath.replace('features.tsv', 'genes.tsv'))
        return compile_mtx(src, prefix, props, clog, raw)

    else:
        return None, [
            f'illegal features.tsv', 
            f'{len(gfile.columns)} columns is observed.'
        ] + [f'[{_x + 1}] ' + str(gfile[_x].tolist()[0:5]) for _x in range(len(gfile.columns))]

    from sh import gzip

    if os.path.exists(f'{src}/{prefix}matrix.mtx'):
        gzip(f'{src}/{prefix}matrix.mtx')
    if os.path.exists(f'{src}/{prefix}features.tsv'):
        gzip(f'{src}/{prefix}features.tsv')
    if os.path.exists(f'{src}/{prefix}barcodes.tsv'):
        gzip(f'{src}/{prefix}barcodes.tsv')
        
    try:
        adata = sc.read_10x_mtx(
            src, var_names = 'gene_ids',
            gex_only = True, make_unique = True,
            prefix = prefix
        )
    except:
        return None, ['read mtxz error.']

    final = process_matrix(props, clog, adata, force_filter = raw)
    
    del adata
    return final, []

def compile_h5(src, h5name, props, clog: textio, raw = False):

    try:
        adata = sc.read_10x_h5(
            f'{src}/{h5name}',
            gex_only = True
        )
        
        if 'gene_ids' not in adata.var.columns:
            del adata
            return None, ['there is no "gene_ids" column in variables metadata.']
    except:
        return None, ['read h5 error']

    lst = adata.var['gene_ids'].tolist()
    lsstab = [y.split('.')[0] for y in lst]
    adata.var_names = lsstab
    del adata.var

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    final = process_matrix(props, clog, adata, force_filter = raw)
    
    del adata
    return final, []
