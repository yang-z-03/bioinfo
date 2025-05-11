
import os
import subprocess
from time import sleep
from printtbl import printtbl
from ansi import common_length
from contextlib import redirect_stdout
import datetime

def status_download(registry):

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
        
        if '.done' in contents: stat += ['done']
        elif '.ignore' in contents: stat += ['ignored']
        elif len(visible) > 0: stat += ['occupied']
        else: stat += ['empty']
    
    return stat

def styles_download(x: str):
    if x == 'done':
        return '\033[32;1m' + x + '\033[0m'
    elif x == 'occupied':
        return '\033[34;1m' + x + '\033[0m'
    elif x == 'empty':
        return '\033[31;1m' + x + '\033[0m'
    else: return '\033[31;0m' + x + '\033[0m'


def status_data(registry):

    # get the download status as a string array
    # which will be appended to the list display.

    stat = []
    data = os.listdir('data')
    for i in range(len(registry['name'])):
        name = registry['name'][i]
        if not name in data:
            stat += ['.']
            continue

        contents = os.listdir('data/{}'.format(name))
        visible = []
        for x in contents:
            if not x.startswith('.'): visible += [x]
        
        if 'data.h5ad' in contents: stat += ['done']
        elif 'src' in contents: stat += ['ready']
        else: stat += ['empty']
    
    return stat

def styles_data(x: str):
    if x == 'done':
        return '\033[32;1m' + x + '\033[0m'
    elif x == 'ready':
        return '\033[34;1m' + x + '\033[0m'
    elif x == 'empty':
        return '\033[31;1m' + x + '\033[0m'
    else: return '\033[31;0m' + x + '\033[0m'


def touch(x):
    if not os.path.exists(x):
        with open(x, 'w') as f:
            os.utime(x, None)
    else: os.utime(x, None)


def display(registry, n_row, show_status = False):

    registry['download'] = status_download(registry)
    registry['data'] = status_data(registry)
    from list import styles_acctype, styles_number
    total_row = len(registry['data'])
    nfrom = 0

    while nfrom < total_row:

        cols = printtbl(
            registry, styles = { 
                'type': styles_acctype, 'accession': styles_number,
                'download': styles_download, 'data': styles_data
            },
            appends = { 
                'type': 9, 'accession': 11,
                'download': 11, 'data': 11
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

        dnstat = registry['download'][nfrom : min(total_row, nfrom + n_row)]
        acctype = registry['type'][nfrom : min(total_row, nfrom + n_row)]
        acc = registry['accession'][nfrom : min(total_row, nfrom + n_row)]
        nm = registry['name'][nfrom : min(total_row, nfrom + n_row)]

        nfrom += n_row
        for dn, acct, acc, name in zip(dnstat, acctype, acc, nm):

            if acct != 'gse':
                print_final('not supported.')
                continue

            if dn == '.':
                os.makedirs('src/{}'.format(name))

            if dn == '.' or dn == 'empty':

                # here, we will retrieve records from geo database, and get the
                # download links for wget to download.

                import GEOparse as geo
                print_replace(common_length('querying geo record ...', 40))
                gse = geo.get_GEO('GSE{}'.format(acc), destdir = 'geo', silent = True)
                
                if 'supplementary_file' not in gse.metadata.keys():
                    continue
                
                files = gse.metadata['supplementary_file']
                downlog = open(f'processed/logs/{name}.dnlog', 'w')
                downlog.writelines([f'job start [{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}] \n\n'])
                downlog.close()
                downlog = open(f'processed/logs/{name}.dnlog', 'a')

                import wget
                for supp in files:

                    fname = supp.split('/')[-1]
                    print_replace(common_length(f'{fname}', 40))

                    retry = 30
                    for xretry in range(retry):

                        def prog(current, total, width = 40):
                            print_replace(
                                common_length(f'{(100 * current / total):.1f}% [{xretry}] ', 10) +
                                common_length(f'{fname}', 40)
                            )

                        try:
                            # try to format the http download link ...
                            httplink = f'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE{acc}&format=file&file={fname}'
                            # wget.download(supp, f'src/{name}/{fname}', bar = prog)
                            wget.download(httplink, f'src/{name}/{fname}', bar = prog)
                            
                            if fname.endswith('_RAW.tar'):
                                print_replace(common_length('extracting ...', 40))
                                import tarfile
                                with tarfile.open(f'src/{name}/{fname}') as tfile:
                                    tarnames = tfile.getnames()
                                    for tn in tarnames:
                                        tfile.extract(tn, f'src/{name}/')
                                        print_replace(common_length(f'extracting {tn} ...', 40))
                                os.remove(f'src/{name}/{fname}')
                            
                            break

                        except:
                            # clean up the debris. and retry download.
                            flist = os.listdir(f'src/{name}')
                            for xf in flist:
                                if xf.startswith(fname + '.'): os.remove(f'src/{name}/{xf}')

                touch(f'src/{name}/.done')
                downlog.close()
                print_final(common_length('done.', 40))
                continue

            print_final(common_length('skipped.', 40))
            continue

        print('\r', end = '\n')
