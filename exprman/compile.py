
import os
import omicverse as ov
import scanpy as sc
import anndata as ad
import pandas as pd
from printtbl import printtbl
from ansi import common_length
from ansi import error, warning, info
from merge import hsize
import datetime
import pickle
from routines import compile_mtx, compile_mtxz, compile_h5

cwd = os.getcwd()
_fe = open('configs/compile.err', 'w')
_fe.writelines([f'[start] {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")} \n\n'])

def display(registry, n_row, show_status = False):

    import warnings

    # ignore all warnings
    warnings.filterwarnings('ignore')
    sc.settings.verbosity = 0

    from download import status_download, styles_download
    from download import status_data, styles_data
    from metagen import status_dtype, styles_dtype

    registry['download'] = status_download(registry)
    registry['dtype'] = status_dtype(registry)
    registry['data'] = status_data(registry)
    from list import styles_acctype, styles_number

    total_row = len(registry['dtype'])
    nfrom = 0

    # here, we will read the sample config file.
    # if the file exist, we will interpret all existing record, and append
    # non-existing datasets to the end of the file. but not changing the previous
    # contents. the sample config is capable for commenting using '#'.
    
    config = pd.read_table('configs/lookup.tsv')
    n_cells = 0
        
    while nfrom < total_row:

        cols = printtbl(
            registry, styles = { 
                'type': styles_acctype, 'accession': styles_number,
                'download': styles_download, 'dtype': styles_dtype,
                'data': styles_data
            },
            appends = { 
                'type': 9, 'accession': 11,
                'download': 11, 'dtype': 11, 'data': 11
            }, nfrom = nfrom, nto = nfrom + n_row - 1
        )

        # move the cursor to the first item.
        from ansi import ansi_move_cursor
        ansi_move_cursor(-(min(total_row, nfrom + n_row) - nfrom), cols)

        def print_final(x):
            print(x, end = '')
            ansi_move_cursor(1, -len(x))

        def print_replace(x):
            print(x, end = '')
            ansi_move_cursor(0, -len(x))
        
        def fprintc(x, l = 40): print_final(common_length(x, l))
        def printc(x, l = 40): print(common_length(x, l), end = '')
        def rprintc(x, l = 40): print_replace(common_length(x, l))

        for s_download, s_name, s_acctype, s_acc, s_data, s_dtype in \
            zip(
                registry['download'][nfrom : min(total_row, nfrom + n_row)], 
                registry['name'][nfrom : min(total_row, nfrom + n_row)], 
                registry['type'][nfrom : min(total_row, nfrom + n_row)], 
                registry['accession'][nfrom : min(total_row, nfrom + n_row)], 
                registry['data'][nfrom : min(total_row, nfrom + n_row)], 
                registry['dtype'][nfrom : min(total_row, nfrom + n_row)]
            ):

            # we read the dimensions, and ensure the indexing to be unique.
            if s_data in ['done']: 
                
                hrsize = hsize(os.path.getsize(f'data/{s_name}/data.h5ad'))
                if os.path.exists(f'data/{s_name}/data.size'):
                    with open(f'data/{s_name}/data.size', 'rb') as dsizef:
                        dsize = pickle.load(dsizef)
                        printc(f'[obs] \033[32;1m{dsize["n_obs"]}\033[0m * [var] \033[31;1m{dsize["n_vars"]}\033[0m {hrsize}', l = 62)
                        ansi_move_cursor(1, -40)
                        n_cells += dsize['n_obs'] if dsize['n_vars'] > 200 else 0
                    continue
                
                h5ad = sc.read_h5ad(f'data/{s_name}/data.h5ad')
                # h5ad.obs['ubc'] = [f'{s_acctype}:{s_acc}' + ':' + str(x + 1) for x in range(h5ad.n_obs)]
                # h5ad.obs_names = [f'{s_acctype}:{s_acc}' + ':' + str(x + 1) for x in range(h5ad.n_obs)]
                # h5ad.write_h5ad(f'data/{s_name}/data.h5ad')
                printc(f'[obs] \033[32;1m{h5ad.n_obs}\033[0m * [var] \033[31;1m{h5ad.n_vars}\033[0m {hrsize}', l = 62)
                with open(f'data/{s_name}/data.size', 'wb') as dsizef:
                    pickle.dump({ 'n_obs': h5ad.n_obs, 'n_vars': h5ad.n_vars }, dsizef)
                ansi_move_cursor(1, -40)
                n_cells += h5ad.n_obs if h5ad.n_vars > 200 else 0
                del h5ad
                continue
                
            if s_dtype in ['unknown', '.', 'empty']: 
                fprintc('data type not supported.')
                continue

            # compilation selection

            contents = os.listdir('src/{}'.format(s_name))
            visible = []
            for x in contents:
                if not x.startswith('.'): visible += [x]

            n_mtx_genes = []
            n_mtx_features = []
            n_h5 = []; f_h5 = []

            for x in visible:
                if x.endswith('matrix.mtx') or x.endswith('matrix.mtx.gz'):
                    root = x.replace('matrix.mtx.gz', '').replace('matrix.mtx', '')
                    if (root + 'genes.tsv') in visible or (root + 'genes.tsv.gz') in visible:
                        n_mtx_genes += [root]
                        continue 
                    elif (root + 'features.tsv') in visible or (root + 'features.tsv.gz') in visible:
                        n_mtx_features += [root]
                        continue
                elif (x.endswith('.h5') or x.endswith('.h5.gz')):
                    f_h5 += [x]
                    n_h5 += [x.replace('.h5.gz', '').replace('.h5', '')]
                
            # create the data directory structures if not exist.
            if s_data == '.':
                os.makedirs(f'data/{s_name}')
                os.symlink(f'{cwd}/src/{s_name}', f'{cwd}/data/{s_name}/src')

            # here, we will set a common interface for the compilation routine
            # and selectively pick one of the routines (based on experience of
            # given data) for analyses.

            # configure properties
                
            props_items = [
                'accession', 'dataset', 'name',
                'tissue', 'genetic', 'strain', 'ms', 'fs', 'sort', 'expm', 'age', 'sex'
            ]

            # filter acceptable samples from the lookup table.
            cfg = config.loc[config['accession'] == 'gse:' + s_acc, :]

            # filter dataset that is single cell rna sequencing but not others.
            is_sc_rna = [((x == 'nan') or ('scrna-gex' in x)) for x in cfg['dstype'].astype('str').tolist()]
            cfg = cfg.loc[is_sc_rna, :]

            # for mus musculus
            cfg = cfg.loc[cfg['tax'] == '10090', :]
            # is selected due to precedency
            cfg = cfg.loc[cfg['selected'] == 'x', :]

            # other fields are not mandatory. let's check if there is any valid
            # sample for the dataset left.
            if len(cfg) == 0:
                fprintc('no valid sample, skipped.')
                continue

            # set up compilation log.

            clog = open(f'processed/logs/compile/{s_name}.clog', 'w')
            clog.writelines([
                f'compilation log \n',
                f'  [dataset] {s_name} ({s_acctype}, {s_acc}) \n',
                f'  [start] {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")} \n\n',
                f'loading dataset components ... \n\n'
            ])

            samples = {} # dictionary of anndata, with obs filled but var left blank.
            
            pd.set_option('display.max_columns', None)
            pd.set_option('display.max_colwidth', 90)
            cfg['name'] = cfg['name'].astype('str')
            print(cfg, file = clog)
            print('', file = clog)
                
            cfgdict = cfg.to_dict('list')
            for c_accession, c_dataset, c_dtype, c_filt, c_encode, c_name, c_prefix, \
                c_strain, c_sex, c_tissue, c_age, c_expm, c_genetic, c_sort in zip(
                    cfgdict['accession'], cfgdict['dataset'], cfgdict['dtype'], cfgdict['filtered'],
                    cfgdict['encode'], cfgdict['name'], cfgdict['prefix'], 
                    cfgdict['strain'], cfgdict['sex'], cfgdict['tissue'], cfgdict['age'],
                    cfgdict['expm'], cfgdict['genetic'], cfgdict['sort']
                ):

                props = {
                    'accession': c_accession,
                    'dataset': c_dataset,
                    'name': c_name if c_name.strip() not in ['nan', ''] else c_encode,
                    'tissue': c_tissue,
                    'genetic': c_genetic,
                    'strain': c_strain,
                    'ms': '.', 'fs': '.',
                    'sort': c_sort,
                    'expm': c_expm,
                    'age': c_age,
                    'sex': c_sex
                }

                raw = (c_filt != "filt")
                sampobj = None
                err = []
                sample = c_prefix if c_prefix not in ['nan'] else ''

                print(f'[{props["accession"]}/{props["name"]}/{c_dtype}]', file = clog)

                try:
                    if c_dtype == "mtx":
                        rprintc(f'{n_cells} [{props["name"]}/mtx]')
                        sampobj, err = compile_mtx(f'src/{s_name}', sample, props, clog, raw = raw)
                    
                    elif c_dtype == "mtxz":
                        rprintc(f'{n_cells} [{props["name"]}/mtxz]')
                        sampobj, err = compile_mtxz(f'src/{s_name}', sample, props, clog, raw = raw)
    
                    elif c_dtype == 'h5':
                        if os.path.exists(f'src/{s_name}/{sample}.h5'):
                            rprintc(f'{n_cells} [{props["name"]}/h5]')
                            sampobj, err = compile_h5(f'src/{s_name}', sample + '.h5', props, clog, raw = raw)
                        elif os.path.exists(f'src/{s_name}/{sample}.h5.gz'):
                            rprintc(f'{n_cells} [{props["name"]}/h5]')
                            sampobj, err = compile_h5(f'src/{s_name}', sample + '.h5.gz', props, clog, raw = raw)
                except:
                    sampleobj = None
                    err = ['compilation or format error breaks the script.']

                if sampobj is not None:
                    samples[props['accession'] + ':' + props['name']] = sampobj
                else: 
                    err_msg = [ f'[{props["accession"]}/{props["name"]}/{c_dtype}] error: \n' ] + \
                        ['    [e] ' + _x + '\n' for _x in err]
                    _fe.writelines(err_msg)
                    clog.writelines(err_msg)
                
                pass

            if len(samples) == 0:
                fprintc('data type totally err.')
                continue
            
            # merge samples. the compile_* methods read files into an data of
            # anndata type, and normalized genes to gene uid, with attached sample
            # metadata onto the obs slot. the var slot is intentionally left blank
            # and left for querying after merging.

            adata = ad.concat(samples, label = "batch", join = 'outer')

            # the merge introduce duplicated names. here, we re-index it.
                
            adata.obs['ubc'] = [f'{s_acctype}:{s_acc}' + ':' + str(x + 1) for x in range(adata.n_obs)]
            adata.obs_names = [f'{s_acctype}:{s_acc}' + ':' + str(x + 1) for x in range(adata.n_obs)]
                
            rprintc('writing to disk ...')
            adata.write_h5ad(f'data/{s_name}/data.h5ad')
            hrsize = hsize(os.path.getsize(f'data/{s_name}/data.h5ad'))
            printc(f'[obs] \033[32;1m{adata.n_obs}\033[0m * [var] \033[31;1m{adata.n_vars}\033[0m {hrsize}', l = 62)
            ansi_move_cursor(1, -40)

            with open(f'data/{s_name}/data.size', 'wb') as dsizef:
                pickle.dump({ 'n_obs': adata.n_obs, 'n_vars': adata.n_vars }, dsizef)
                n_cells += adata.n_obs if adata.n_vars > 200 else 0
            
            clog.writelines([
                '\n', 'export data: \n\n', 
                f'  anndata with obs * var: {adata.n_obs} * {adata.n_vars} \n'])
            clog.writelines(['\n', 'observation metadata: \n\n'])
            print(adata.obs, file = clog)

            clog.flush()
            _fe.flush()
            clog.close()

            keys = list(samples.keys())
            for k in keys: del samples[k]
            del samples
            del adata
            os.symlink(f'{cwd}/data/{s_name}/data.h5ad', f'{cwd}/adata/{s_name}.h5ad')
            pass
        
        # update iteration
        nfrom += n_row
        print('\r', end = '\n')
        pass
