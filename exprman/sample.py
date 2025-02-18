
# produce a sample table with linked files.
# we will not use the printtbl routines, but rather print by ourself.

def formal_length(*_x, l = 40):
    
    x = ''
    for _item in _x:
        x += _item
    
    real_len = 0
    in_033_block = False
    
    for c in x:
        if c == '\033':
            in_033_block = True
            continue
        if c == 'm' and in_033_block:
            in_033_block = False
            continue

        real_len += 1
        continue
    
    construct = ''
    if real_len > l:
        wordlen = 0
        for c in x:
            if c == '\033':
                in_033_block = True
                construct += '\033'
                continue
            if c == 'm' and in_033_block:
                in_033_block = False
                construct += 'm'
                continue
            if in_033_block:
                construct += c
                continue
            
            wordlen += 1
            construct += c
            if wordlen >= l - 4:
                if in_033_block: construct += '\033[0m'
                construct += ' ...'
                break
    
    else: construct = x + (l - real_len) * ' '
    return construct

def red(x): return '\033[31m' + x + '\033[0m'
def green(x): return '\033[32m' + x + '\033[0m'
def yellow(x): return '\033[33m' + x + '\033[0m'
def cyan(x): return '\033[36m' + x + '\033[0m'
def printe(*a): print(*a, end = '', sep = '')

def display(registry, n_row, show_status = False):
    
    for s_name, s_type, s_acc in zip(
        registry['name'], 
        registry['type'], 
        registry['accession']
    ):
        
        if s_type != 'gse':
            printe(formal_length(red(s_type), ':', yellow(s_acc), l = 40))
            printe(formal_length('dataset not supported.'))
            print()
            continue

        gserec = 'GSE' + s_acc
        import GEOparse as geo
        gse = geo.get_GEO(gserec, destdir = 'geo', silent = True)

        fcatlog = {
            'dataset': {},
            'samples': {}
        }

        # dataset information

        if 'title' in gse.metadata.keys():
            fcatlog['dataset']['title'] = gse.metadata['title']
        else: fcatlog['dataset']['title'] = None

        if 'geo_accession' in gse.metadata.keys():
            fcatlog['dataset']['acc'] = gse.metadata['geo_accession']
        else: fcatlog['dataset']['acc'] = None

        if 'sample_taxid' in gse.metadata.keys():
            fcatlog['dataset']['tax'] = gse.metadata['sample_taxid']
        else: fcatlog['dataset']['tax'] = None
        
        if 'last_update_date' in gse.metadata.keys():
            fcatlog['dataset']['date'] = gse.metadata['last_update_date']
        else: fcatlog['dataset']['date'] = None

        if 'summary' in gse.metadata.keys():
            fcatlog['dataset']['summary'] = gse.metadata['summary']
        else: fcatlog['dataset']['summary'] = None

        if 'overall_design' in gse.metadata.keys():
            fcatlog['dataset']['design'] = gse.metadata['overall_design']
        else: fcatlog['dataset']['design'] = None

        fcatlog['dataset']['files'] = []
        if 'supplementary_file' in gse.metadata.keys():
            for supp in gse.metadata['supplementary_file']:
                fname = supp.split('/')[-1]
                if fname == gserec + '_RAW.tar': continue
                else: fcatlog['dataset']['files'] += [fname]
                
        # sample information
        
        for k in gse.gsms.keys():
            samp = gse.gsms[k]
            samp_item = {}

            if 'source_name_ch1' in samp.metadata.keys():
                samp_item['name'] = samp.metadata['source_name_ch1']
            else: samp_item['name'] = None

            if 'taxid_ch1' in samp.metadata.keys():
                samp_item['tax'] = samp.metadata['taxid_ch1']
            else: samp_item['tax'] = None
            
            if 'geo_accession' in samp.metadata.keys():
                samp_item['acc'] = samp.metadata['geo_accession']
            else: samp_item['acc'] = None

            if 'last_update_date' in samp.metadata.keys():
                samp_item['date'] = samp.metadata['last_update_date']
            else: samp_item['date'] = None

            if 'characteristics_ch1' in samp.metadata.keys():
                samp_item['char'] = samp.metadata['characteristics_ch1']
            else: samp_item['char'] = None

            if 'treatment_protocol_ch1' in samp.metadata.keys():
                samp_item['prtcl.treat'] = samp.metadata['treatment_protocol_ch1']
            else: samp_item['prtcl.treat'] = None

            if 'extract_protocol_ch1' in samp.metadata.keys():
                samp_item['prtcl.extract'] = samp.metadata['extract_protocol_ch1']
            else: samp_item['prtcl.extract'] = None

            samp_item['files'] = []
            for metakey in samp.metadata.keys():
                if metakey.startswith('supplementary_file'):
                    for supp in samp.metadata[metakey]:
                        fname = supp.split('/')[-1]
                        samp_item['files'] += [fname]
            
            fcatlog['samples'][k] = samp_item

        # list downloaded files
        import os
        down = os.listdir(f'src/{s_name}')
        visible = []
        for x in down:
            if not x.startswith('.'): visible += [x]
        
        n_mtx_genes = []
        n_mtx_features = []
        n_raw_h5 = []
        n_filtered_h5 = []
        n_undef_h5 = []
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
            elif (x.endswith('.h5') or x.endswith('.h5.gz')) and '_raw_' in x:
                n_raw_h5 += [x.replace('raw_feature_bc_matrix.h5.gz', '').replace('raw_feature_bc_matrix.h5', '')]
            elif (x.endswith('.h5') or x.endswith('.h5.gz')) and '_filtered_' in x:
                n_filtered_h5 += [x.replace('filtered_feature_bc_matrix.h5.gz', '').replace('filtered_feature_bc_matrix.h5', '')]
            elif (x.endswith('.h5') or x.endswith('.h5.gz')):
                n_undef_h5 += [x.replace('.h5.gz', '').replace('.h5', '')]
            else: unrecog += [x]
        
        for x in removed:
            unrecog.remove(x)

        # precedency rules
        # h5f > h5r > mtx + mtxz > h5u

        # if sample contains its own set of matrix, we will ignore those provided
        # as a whole in the dataset root. if all samples provide no valid matrix
        # we will search the root directory.

        printe(formal_length(red(s_type), ':', yellow(s_acc), l = 40))
        if len(n_filtered_h5) > 0:
            printe('[h5f] ', red(len(n_filtered_h5)), ' ')
        elif len(n_raw_h5) > 0:
            printe('[h5r] ', red(len(n_raw_h5)), ' ')
        elif len(n_mtx_features) + len(n_mtx_genes) > 0:
            if len(n_mtx_features) > 0:
                printe('[mtxz] ', red(len(n_mtx_features)), ' ')
            if len(n_mtx_genes) > 0:
                printe('[mtx] ', red(len(n_mtx_genes)), ' ')
        elif len(n_undef_h5) > 0:
            printe('[h5u] ', red(len(n_undef_h5)), ' ')
        
        print()

        # find the name from the catlog tree.
        def findname(root, catlog):
            for samplek in catlog['samples'].keys():
                sample = catlog['samples'][samplek]
                for fname in sample['files']:
                    if fname.startswith(root):
                        return {
                            'acc': sample['acc'],
                            'name': sample['name']
                        }
            
            for fname in catlog['dataset']['files']:
                if fname.startswith(root):
                    snames = list(catlog['samples'].keys())
                    for sname in snames:
                        if sname.lower() in fname:
                            return {
                                'acc': catlog['samples'][sname]['acc'],
                                'name': catlog['samples'][sname]['name']
                            }
                    return {
                        'acc': catlog['dataset']['acc'],
                        'name': 'root'
                    }
            return {
                'acc': 'fail',
                'name': 'fail'
            }
        
        def uniform(x):
            return x.lower().replace('-','.').replace('_','.').replace('   ', '.').replace('  ', '.').replace(' ', '.')
        
        def printinfo(prefix, catlog):
            fn = findname(prefix, catlog)
                
            if fn['acc'] == 'fail':
                printe(red('cannot find parent.'))
                return
            elif fn['acc'].lower().startswith('gse'):
                printe(formal_length(red('gse'), ':', yellow(fn['acc'].lower().replace('gse', '')), l = 16))
            elif fn['acc'].lower().startswith('gsm'):
                printe(formal_length(red('gsm'), ':', green(fn['acc'].lower().replace('gsm', '')), l = 16))
            else: printe(formal_length(red('?'), l = 16))

            # print name only
            if fn['acc'] == catlog['dataset']['acc']:
                printe(formal_length('root', l = 30))
            elif fn['acc'] in catlog['samples'].keys():
                printe(formal_length(uniform(catlog['samples'][fn['acc']]['name']), l = 30))
        
        if len(n_filtered_h5) > 0:
            for fp in n_filtered_h5:
                printe(formal_length('| ' + fp))
                printinfo(fp + 'filtered_feature_bc_matrix.h5', fcatlog)
                print()

        elif len(n_raw_h5) > 0:
            for fp in n_raw_h5:
                printe(formal_length('| ' + fp))
                printinfo(fp + 'raw_feature_bc_matrix.h5', fcatlog)
                print()

        elif len(n_mtx_features) + len(n_mtx_genes) > 0:
            if len(n_mtx_features) > 0:
                for fp in n_mtx_features:
                    printe(formal_length('| ' + fp))
                    printinfo(fp + 'matrix.mtx', fcatlog)
                    print()

            if len(n_mtx_genes) > 0:
                for fp in n_mtx_genes:
                    printe(formal_length('| ' + fp))
                    printinfo(fp + 'matrix.mtx', fcatlog)
                    print()

        elif len(n_undef_h5) > 0:
            for fp in n_raw_h5:
                printe(formal_length('| ' + fp))
                printinfo(fp + '.h5', fcatlog)
                print()
