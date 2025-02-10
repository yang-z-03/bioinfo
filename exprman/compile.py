
import os
import omicverse as ov
import anndata as ad
import pandas as pd
from printtbl import printtbl
from ansi import common_length
from ansi import error, warning, info

cwd = os.getcwd()

def display(registry, n_row, show_status = False):

    from download import status_download, styles_download
    from download import status_data, styles_data
    from metagen import status_dtype, styles_dtype

    registry['download'] = status_download(registry)
    registry['dtype'] = status_dtype(registry)
    registry['data'] = status_data(registry)
    from list import styles_acctype, styles_number

    total_row = len(registry['dtype'])
    nfrom = 0

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
        def rprintc(x, l = 40): print_replace(common_length(x, l))

        # here, we will read the sample config file.
        # if the file exist, we will interpret all existing record, and append
        # non-existing datasets to the end of the file. but not changing the previous
        # contents. the sample config is capable for commenting using '#'.

        config = {}

        current_dataset = '*'
        current_sample = '*'

        if os.path.exists('configs/samples'):
            with open('configs/samples', 'r') as fr:
                contents = fr.read().splitlines()
                for line in contents:
                    if len(line.strip()) == 0: continue
                    if line.strip().startswith('#'): continue

                    tokens = line.strip().split()
                    if len(tokens) < 2: 
                        error('syntax error:', line)
                        continue

                    if tokens[0] == 'dataset':
                        current_dataset = tokens[1]
                        current_sample = '*'
                        continue

                    elif tokens[0] == 'sample':
                        current_sample = tokens[1]
                        continue

                    elif tokens[0] == 'prop':
                        if len(tokens) <= 2: continue
                        if current_dataset not in config.keys():
                            config[current_dataset] = {}
                        if current_sample not in config[current_dataset].keys():
                            config[current_dataset][current_sample] = {}
                        config[current_dataset][current_sample][tokens[1]] = tokens[2]

        for s_download, s_name, s_acctype, s_acc, s_data, s_dtype in \
            zip(
                registry['download'][nfrom : min(total_row, nfrom + n_row)], 
                registry['name'][nfrom : min(total_row, nfrom + n_row)], 
                registry['type'][nfrom : min(total_row, nfrom + n_row)], 
                registry['accession'][nfrom : min(total_row, nfrom + n_row)], 
                registry['data'][nfrom : min(total_row, nfrom + n_row)], 
                registry['dtype'][nfrom : min(total_row, nfrom + n_row)]
            ):

            if s_dtype in ['unknown', '.', 'empty']: 
                fprintc('data type not supported.')
                continue
            if s_data in ['done']: 
                fprintc('compilation done, skipped.')
                continue

            # compilation selection

            contents = os.listdir('src/{}'.format(s_name))
            visible = []
            for x in contents:
                if not x.startswith('.'): visible += [x]

            n_mtx_genes = []
            n_mtx_features = []
            n_raw_h5 = []; f_raw_h5 = []
            n_filtered_h5 = []; f_filtered_h5 = []
            n_undef_h5 = []; f_undef_h5 = []

            for x in visible:
                if x.endswith('matrix.mtx') or x.endswith('matrix.mtx.gz'):
                    root = x.replace('matrix.mtx.gz', '').replace('matrix.mtx', '')
                    if (root + 'genes.tsv') in visible or (root + 'genes.tsv.gz') in visible:
                        n_mtx_genes += [root]
                        continue 
                    elif (root + 'features.tsv') in visible or (root + 'features.tsv.gz') in visible:
                        n_mtx_features += [root]
                        continue
                elif (x.endswith('.h5') or x.endswith('.h5.gz')) and '_raw_' in x:
                    f_raw_h5 += [x]
                    n_raw_h5 += [x.replace('raw_feature_bc_matrix.h5.gz', '').replace('raw_feature_bc_matrix.h5', '')]
                elif (x.endswith('.h5') or x.endswith('.h5.gz')) and '_filtered_' in x:
                    f_filtered_h5 += [x]
                    n_filtered_h5 += [x.replace('filtered_feature_bc_matrix.h5.gz', '').replace('filtered_feature_bc_matrix.h5', '')]
                elif (x.endswith('.h5') or x.endswith('.h5.gz')):
                    f_undef_h5 += [x]
                    n_undef_h5 += [x.replace('.h5.gz', '').replace('.h5', '')]
                
            # create the data directory structures if not exist.
            if s_data == '.':
                os.makedirs(f'data/{s_name}')
                os.symlink(f'{cwd}/src/{s_name}', f'{cwd}/data/{s_name}/src')

            # here, we will set a common interface for the compilation routine
            # and selectively pick one of the routines (based on experience of
            # given data) for analyses.

            # configure properties

            props_items = ['name', 'tissue', 'genetic', 'strain', 'sort', 'expm', 'age', 'sex']
            def get_props(s_sample):
                
                pdict = {}
                error_code = 0
                for x in props_items:
                    if s_name not in config.keys():
                        error_code = 1
                        break
                    if s_sample not in config[s_name].keys():
                        error_code = 2
                        break

                    if x in config[s_name][s_sample].keys():
                        pdict[x] = config[s_name][s_sample][x]
                        continue
                    
                    if '*' in config[s_name].keys():
                        if x in config[s_name]['*'].keys():
                            pdict[x] = config[s_name]['*'][x]
                            continue
                    
                    if '*' in config.keys():
                        if '*' in config['*'].keys():
                            if x in config['*']['*'].keys():
                                pdict[x] = config['*']['*'][x]
                                continue
                    
                    error_code = 3
                    break
                
                pdict['accession'] = s_acctype + ':' + str(s_acc)
                pdict['dataset'] = s_name.replace(str(s_acc) + '-', '').replace(str(s_acc), '')
                pdict['ms'] = '.'
                pdict['fs'] = '.'

                return pdict, error_code
            
            if s_name not in config.keys():
                
                # do not find a valid configuration.
                # just skip it without any consideration.

                fprintc('no valid config, skipped.')
                continue

            # set up compilation log.

            import datetime

            clog = open(f'processed/logs/{s_name}.clog', 'w')
            clog.writelines([
                f'compilation log \n',
                f'  [dataset] {s_name} ({s_acctype}, {s_acc}) \n',
                f'  [start] {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")} \n\n',
                f'loading dataset components ... \n'
            ])

            from routines import compile_mtx, compile_mtxz, compile_h5r, compile_h5f
            samples = {} # dictionary of anndata, with obs filled but var left blank.

            if len(n_mtx_genes) + len(n_mtx_features) > 0:

                for sample in n_mtx_genes: 
                    props, e = get_props(sample)
                    if e != 0: 
                        clog.writelines([
                            f'  [e] sample [{sample}/mtx] skipped due to illegal property set, error code = {e}. \n'
                        ])
                        continue
                    clog.writelines(
                        [f'  [i] sample: {sample} \n'] + 
                        [f'    [i] property {z}: {props[z]} \n' for z in props.keys()]
                    )
                    rprintc(f'reading [{props["name"]}/mtx]')
                    samples[props['accession'] + ':' + props['name']] = \
                        compile_mtx(f'src/{s_name}', sample, props, clog)

                for sample in n_mtx_features: 
                    props, e = get_props(sample)
                    if e != 0: 
                        clog.writelines([
                            f'  [e] sample [{sample}/mtxz] skipped due to illegal property set, error code = {e}. \n'
                        ])
                        continue
                    clog.writelines(
                        [f'  [i] sample: {sample} \n'] + 
                        [f'    [i] property {z}: {props[z]} \n' for z in props.keys()]
                    )
                    rprintc(f'reading [{props["name"]}/mtxz]')
                    samples[props['accession'] + ':' + props['name']] = \
                        compile_mtxz(f'src/{s_name}', sample, props, clog)
            
            elif len(n_filtered_h5) > 0:

                for sample, samplef in zip(n_filtered_h5, f_filtered_h5): 
                    props, e = get_props(sample)
                    if e != 0: 
                        clog.writelines([
                            f'  [e] sample [{sample}/h5f] skipped due to illegal property set, error code = {e}. \n'
                        ])
                        continue
                    clog.writelines(
                        [f'  [i] sample: {sample} \n'] + 
                        [f'    [i] property {z}: {props[z]} \n' for z in props.keys()]
                    )
                    rprintc(f'reading [{props["name"]}/h5f]')
                    samples[props['accession'] + ':' + props['name']] = \
                        compile_h5f(f'src/{s_name}', samplef, props, clog)
            
            elif len(n_raw_h5) + len(n_undef_h5) > 0:
                
                for sample, samplef in zip(n_raw_h5, f_raw_h5): 
                    props, e = get_props(sample)
                    if e != 0: 
                        clog.writelines([
                            f'  [e] sample [{sample}/h5ru] skipped due to illegal property set, error code = {e}. \n'
                        ])
                        continue
                    clog.writelines(
                        [f'  [i] sample: {sample} \n'] + 
                        [f'    [i] property {z}: {props[z]} \n' for z in props.keys()]
                    )
                    rprintc(f'reading [{props["name"]}/h5ru]')
                    samples[props['accession'] + ':' + props['name']] = \
                        compile_h5r(f'src/{s_name}', samplef, props, clog)

                for sample, samplef in zip(n_undef_h5, f_undef_h5): 
                    props, e = get_props(sample)
                    if e != 0: 
                        clog.writelines([
                            f'  [e] sample [{sample}/h5ru] skipped due to illegal property set, error code = {e}. \n'
                        ])
                        continue
                    clog.writelines(
                        [f'  [i] sample: {sample} \n'] + 
                        [f'    [i] property {z}: {props[z]} \n' for z in props.keys()]
                    )
                    rprintc(f'reading [{props["name"]}/h5ru]')
                    samples[props['accession'] + ':' + props['name']] = \
                        compile_h5r(f'src/{s_name}', samplef, props, clog)

            # merge samples. the compile_* methods read files into an data of
            # anndata type, and normalized genes to gene uid, with attached sample
            # metadata onto the obs slot. the var slot is intentionally left blank
            # and left for querying after merging.

            adata = ad.concat(samples, label = "batch", join = 'outer')
            rprintc('writing to disk ...')
            adata.write_h5ad(f'data/{s_name}/data.h5ad')
            fprintc('done.')

            pd.set_option('display.max_columns', None)
            pd.set_option('display.max_colwidth', 90)
            clog.writelines([
                '\n', 'export data: \n\n', 
                f'  anndata with obs * var: {adata.n_obs} * {adata.n_vars} \n'])
            clog.writelines(['\n', 'observation metadata: \n\n'])
            print(adata.obs, file = clog)

            clog.flush()
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
