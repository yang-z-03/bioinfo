
import curses
import argparse
import os

from ansi import error, warning, info

parser = argparse.ArgumentParser(
    prog = 'exprman',
    description = 'manage expression matrix and dataset integration'
)

parser.add_argument('dir', type = str, help = 'working directory')
parser.add_argument('-n', dest = 'nrow', default = 10, type = int, help = 'number of rows per display')
parser.add_argument('-t', dest = 'task', default = 'list', type = str, help = 'display task for manager')

# clean options
parser.add_argument('--clean', dest = 'clean', default = '', type = str, help = 'which partition to clean up')
parser.add_argument('--accession', dest = 'acc', default = '', type = str, help = 'dataset accession to clean up')
parser.add_argument('--from', dest = 'ifrom', default = 9999, type = int, help = 'dataset accession from index to clean up')
parser.add_argument('--to', dest = 'ito', default = 9999, type = int, help = 'dataset accession to index to clean up')
parser.add_argument('--dry', dest = 'dry', action = 'store_true', default = False, help = 'dry run the cleaning process')

parser.add_argument('--status', action = 'store_true', default = False, help = 'print status message')

opt = parser.parse_args()
try: os.chdir(opt.dir)
except: error(f'{opt.dir} not redirectable.')

# test if the required dataset catalog exists, if exist, set up the sample table
# automatically for non-configured projects.

if not os.path.exists('configs/register'):
    error('could not find configs/register.')

register = []
with open('configs/register', 'r') as f:
    register = f.read().splitlines()

reg = {
    'type': [],
    'accession': [],
    'name': []
}

elem0s = []
print('')
has_warn = 0
for line in register:
    if line.strip() == '': continue
    if line.strip().startswith('#'): continue
    
    elems = line.split()

    if elems[0] in elem0s:
        warning(f'duplicated accession [{elems[0]}]. discarding name [{elems[1]}].')
        has_warn += 1
        continue

    reg['name'] += [elems[1]]
    _type, _acc = elems[0].split(':')
    reg['type'] += [_type]
    reg['accession'] += [_acc]
    elem0s += [elems[0]]

if has_warn > 0: print('')

# filter dataset.
reg['type'] = reg['type'][opt.ifrom:opt.ito]
reg['accession'] = reg['accession'][opt.ifrom:opt.ito]
reg['name'] = reg['name'][opt.ifrom:opt.ito]

# length of registered datasets.
lenreg = len(reg['name'])

if opt.task == 'list':
    from list import display
    display(reg, opt.nrow, opt.status)

elif opt.task == 'download':
    from download import display
    display(reg, opt.nrow, opt.status)

elif opt.task == 'metagen':
    from metagen import display
    display(reg, opt.nrow, opt.status)

elif opt.task == 'metachk':
    from metachk import display
    display(reg, opt.nrow, opt.status)

elif opt.task == 'compile':
    from compile import display
    display(reg, opt.nrow, opt.status)

elif opt.task == 'qc':
    from qc import display
    display(reg, opt.nrow, opt.status)

elif opt.task == 'merge':
    from merge import display
    display(reg, opt.nrow, opt.status)

elif opt.task == 'sample':
    from sample import display
    display(reg, opt.nrow, opt.status)

elif opt.task == 'chat':
    import chat

elif opt.task == 'clean':
    searchs = [x + ':' + y for x, y in zip(reg['type'], reg['accession'])]
    names = reg['name']
    
    dsnames = names
    
    if opt.acc in searchs: 
        dsnames = [names[searchs.index(opt.acc)]]
        info(f'cleaning up [{opt.clean}] partion for dataset [{opt.acc}] {dsnames[0]}')
    else: pass

    if len(dsnames) == 0:
        error('unrecognized accession number. check configs/register for details')
        
    for dsname in dsnames:
        
        if opt.clean in ['qc', 'compile', 'download']:

            import shutil
            for fname in os.listdir('processed/graphs/detections'):
                if fname == dsname.lower():
                    if not opt.dry: shutil.rmtree('processed/graphs/detections/' + fname.lower())
                    warning(f'removed processed/graphs/detections/{fname.lower()}/*')
            
            for fname in os.listdir('processed/graphs/mito'):
                if fname == dsname.lower():
                    if not opt.dry: shutil.rmtree('processed/graphs/mito/' + fname.lower())
                    warning(f'removed processed/graphs/mito/{fname.lower()}/*')
            
            for fname in os.listdir('processed/hvg'):
                if fname.startswith(dsname + '.'):
                    if not opt.dry: os.remove('processed/hvg/' + fname)
                    warning(f'removed processed/hvg/{fname}')
            
            for fname in os.listdir('processed/logs/qc'):
                if fname.startswith(dsname + '.qclog'):
                    if not opt.dry: os.remove('processed/logs/qc/' + fname)
                    warning(f'removed processed/logs/{fname}')
            
            for fname in os.listdir('processed/metrics'):
                if fname.startswith(dsname + '.'):
                    if not opt.dry: os.remove('processed/metrics/' + fname)
                    warning(f'removed processed/metrics/{fname}')
            
            for fname in os.listdir('processed/qc'):
                if fname.startswith(dsname + '.'):
                    if not opt.dry: os.remove('processed/qc/' + fname)
                    warning(f'removed processed/qc/{fname}')
        
        if opt.clean in ['compile', 'download']:
    
            for fname in os.listdir('adata'):
                if fname.startswith(dsname + '.'):
                    if not opt.dry: os.remove('adata/' + fname)
                    warning(f'unlinked adata/{fname}')
            
            for fname in os.listdir('data'):
                if fname == dsname:
                    if os.path.exists(f'data/{fname}/data.h5ad'):
                        if not opt.dry: os.remove(f'data/{fname}/data.h5ad')
                        warning(f'removed data/{fname}/data.h5ad')
        
        if opt.clean in ['download']:
            
            for fname in os.listdir('data'):
                if fname.startswith(dsname):
                    if os.path.exists(f'data/{fname}/src'):
                        if not opt.dry: os.remove(f'data/{fname}/src')
                        warning(f'unlinked data/{fname}/src')
            
            for fname in os.listdir('src'):
                if fname == dsname:
                    for f in os.listdir('src/' + fname):
                        if not opt.dry: os.remove(f'src/{fname}/{f}')
                        warning(f'removed src/{fname}/{f}')

print('')
