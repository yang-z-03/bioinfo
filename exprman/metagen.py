
import os
import subprocess
from time import sleep
from printtbl import printtbl
from ansi import common_length
from ansi import error

def status_dtype(registry):

    # get the download status as a string array
    # which will be appended to the list display.

    stat = []
    src = os.listdir('src')
    for i in range(len(registry['name'])):
        name = registry['name'][i]
        if not name in src:
            stat += ['.']
            continue

        contents = os.listdir('src/{}'.format(name))
        visible = []
        for x in contents:
            if not x.startswith('.'): visible += [x]
        
        n_mtx_genes = 0
        n_mtx_features = 0
        n_h5 = 0

        for x in visible:
            if x.endswith('matrix.mtx') or x.endswith('matrix.mtx.gz'):
                root = x.replace('matrix.mtx.gz', '').replace('matrix.mtx', '')
                if (root + 'genes.tsv') in visible or (root + 'genes.tsv.gz') in visible:
                    n_mtx_genes += 1
                    continue 
                elif (root + 'features.tsv') in visible or (root + 'features.tsv.gz') in visible:
                    n_mtx_features += 1
                    continue
            elif (x.endswith('.h5') or x.endswith('.h5.gz')):
                n_h5 += 1
        
        qstr = f'mtx:{n_mtx_genes} ' if n_mtx_genes > 0 else ''
        qstr += f'mtxz:{n_mtx_features} ' if n_mtx_features > 0 else ''
        qstr += f'h5:{n_h5}' if n_h5 > 0 else ''
        qstr = qstr.strip()

        if n_mtx_genes + n_mtx_features + n_h5 > 0: stat += [qstr]
        elif len(visible) > 0: stat += ['unknown']
        else: stat += ['empty']
    
    return stat

def styles_dtype(x: str):
    if x == 'unknown':
        return '\033[31;1m' + x + '\033[0m'
    elif x == 'empty' or x == '.':
        return '\033[34;1m' + x + '\033[0m'
    else: return '\033[31;0m' + x + '\033[0m'


def touch(x):
    if not os.path.exists(x):
        with open(x, 'w') as f:
            os.utime(x, None)
    else: os.utime(x, None)


def display(registry, n_row, show_status = False):

    from download import status_download, styles_download
    from download import status_data

    registry['download'] = status_download(registry)
    registry['dtype'] = status_dtype(registry)
    from list import styles_acctype, styles_number

    cols = printtbl(
        registry, styles = { 
            'type': styles_acctype, 'accession': styles_number,
            'download': styles_download, 'dtype': styles_dtype
        },
        appends = { 
            'type': 9, 'accession': 11,
            'download': 11, 'dtype': 11
        }
    )

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

    # find missing datasets
    added_lines = []
    registry['data'] = status_data(registry)

    for s_download, s_name, s_acctype, s_acc, s_data, s_dtype in \
        zip(
            registry['download'], 
            registry['name'], 
            registry['type'], 
            registry['accession'], 
            registry['data'], 
            registry['dtype']
        ):
        
        if s_dtype in ['unknown', '.', 'empty']: continue
        if s_data in ['done']: continue

        contents = os.listdir('src/{}'.format(s_name))
        visible = []
        for x in contents:
            if not x.startswith('.'): visible += [x]
        
        n_mtx_genes = []
        n_mtx_features = []
        n_h5 = []
        unrecog = []
        removed = []

        for x in visible:
            if x.endswith('matrix.mtx') or x.endswith('matrix.mtx.gz'):
                root = x.replace('matrix.mtx.gz', '').replace('matrix.mtx', '')
                if (root + 'genes.tsv') in visible or (root + 'genes.tsv.gz') in visible:
                    n_mtx_genes += [root]
                    if (root + 'genes.tsv') in visible: removed.append(root + 'genes.tsv')
                    if (root + 'genes.tsv.gz') in visible: removed.append(root + 'genes.tsv.gz')
                    if (root + 'barcodes.tsv') in visible: removed.append(root + 'barcodes.tsv')
                    if (root + 'barcodes.tsv.gz') in visible: removed.append(root + 'barcodes.tsv.gz')
                    continue 
                elif (root + 'features.tsv') in visible or (root + 'features.tsv.gz') in visible:
                    n_mtx_features += [root]
                    if (root + 'features.tsv') in visible: removed.append(root + 'features.tsv')
                    if (root + 'features.tsv.gz') in visible: removed.append(root + 'features.tsv.gz')
                    if (root + 'barcodes.tsv') in visible: removed.append(root + 'barcodes.tsv')
                    if (root + 'barcodes.tsv.gz') in visible: removed.append(root + 'barcodes.tsv.gz')
                    continue
            elif (x.endswith('.h5') or x.endswith('.h5.gz')):
                n_h5 += [x.replace('.h5.gz', '').replace('.h5', '')]
            else: unrecog += [x]
        
        for x in removed:
            unrecog.remove(x)
        
        default_props = [
            f'        prop name    ',
            f'        prop tissue  ',
            f'        prop genetic ', 
            f'        prop strain  ',
            f'        prop sort    ',
            f'        prop expm    ',
            f'        prop age     ',
            f'        prop sex     '
        ]

        # no previous config.
        if s_name not in config.keys():

            added_lines += [f'dataset {s_name} {s_acctype}:{s_acc}']
            for samp in n_mtx_genes:
                added_lines += [f'    sample {samp} [mtx]'] + default_props
            for samp in n_mtx_features:
                added_lines += [f'    sample {samp} [mtxz]'] + default_props
            for samp in n_h5:
                added_lines += [f'    sample {samp} [h5u]'] + default_props
            
            for files in unrecog:
                added_lines += [f'    unrecog {files}']
            added_lines += ['']
        
    with open('configs/samples', 'a') as fa:
        fa.writelines([x + '\n' for x in added_lines])
    
    pass
