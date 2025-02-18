
# produce a sample table with linked files.
# we will not use the printtbl routines, but rather print by ourself.

_f = open('configs/lookup', 'w')

def formal_length(*_x, l = 40):
    
    x = ''
    for _item in _x:
        x += str(_item)
    
    real_len = 0
    in_033_block = False
    
    for c in x:
        if c == '\033':
            in_033_block = True
            continue
        if c == 'm' and in_033_block:
            in_033_block = False
            continue
        if in_033_block:
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
    
    else: construct = x + ((l - real_len) * ' ')
    return construct

def red(x): return '\033[31m' + str(x) + '\033[0m'
def green(x): return '\033[32m' + str(x) + '\033[0m'
def yellow(x): return '\033[33m' + str(x) + '\033[0m'
def cyan(x): return '\033[36m' + str(x) + '\033[0m'
def printe(*a): print(*a, end = '', sep = '')

def display(registry, n_row, show_status = False):

    n_total = 0
    n_available_ds = 0
    n_available_sample = 0
    
    for s_name, s_type, s_acc in zip(
        registry['name'], 
        registry['type'], 
        registry['accession']
    ):

        n_total += 1
        if s_type != 'gse':
            printe(formal_length(red(s_type), ':', yellow(s_acc), l = 40))
            printe(formal_length(' dataset not supported.'))
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
        else: fcatlog['dataset']['title'] = [None]

        if 'geo_accession' in gse.metadata.keys():
            fcatlog['dataset']['acc'] = gse.metadata['geo_accession']
        else: fcatlog['dataset']['acc'] = [None]

        if 'sample_taxid' in gse.metadata.keys():
            fcatlog['dataset']['tax'] = gse.metadata['sample_taxid']
        else: fcatlog['dataset']['tax'] = [None]
        
        if 'last_update_date' in gse.metadata.keys():
            fcatlog['dataset']['date'] = gse.metadata['last_update_date']
        else: fcatlog['dataset']['date'] = [None]

        if 'summary' in gse.metadata.keys():
            fcatlog['dataset']['summary'] = gse.metadata['summary']
        else: fcatlog['dataset']['summary'] = []

        if 'overall_design' in gse.metadata.keys():
            fcatlog['dataset']['design'] = gse.metadata['overall_design']
        else: fcatlog['dataset']['design'] = []

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
            else: samp_item['name'] = [None]
            
            if 'geo_accession' in samp.metadata.keys():
                samp_item['acc'] = samp.metadata['geo_accession']
            else: samp_item['acc'] = [None]

            if 'taxid_ch1' in samp.metadata.keys():
                samp_item['tax'] = samp.metadata['taxid_ch1']
            else: samp_item['tax'] = [None]

            if 'last_update_date' in samp.metadata.keys():
                samp_item['date'] = samp.metadata['last_update_date']
            else: samp_item['date'] = [None]

            if 'characteristics_ch1' in samp.metadata.keys():
                samp_item['char'] = samp.metadata['characteristics_ch1']
            else: samp_item['char'] = []

            if 'treatment_protocol_ch1' in samp.metadata.keys():
                samp_item['prtcl.treat'] = samp.metadata['treatment_protocol_ch1']
            else: samp_item['prtcl.treat'] = []

            if 'extract_protocol_ch1' in samp.metadata.keys():
                samp_item['prtcl.extract'] = samp.metadata['extract_protocol_ch1']
            else: samp_item['prtcl.extract'] = []

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
            elif (x.endswith('.h5') or x.endswith('.h5.gz')) and 'raw_' in x:
                n_raw_h5 += [x.replace('raw_feature_bc_matrix.h5.gz', '').replace('raw_feature_bc_matrix.h5', '') \
                              .replace('raw_gene_bc_matrices_h5.h5.gz', '').replace('raw_gene_bc_matrices_h5.h5', '')]
            elif (x.endswith('.h5') or x.endswith('.h5.gz')) and 'filtered_' in x:
                n_filtered_h5 += [x.replace('filtered_feature_bc_matrix.h5.gz', '').replace('filtered_feature_bc_matrix.h5', '') \
                                   .replace('filtered_gene_bc_matrices_h5.h5.gz', '').replace('filtered_gene_bc_matrices_h5.h5', '')]
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
            printe(' [h5f] ', red(len(n_filtered_h5)), '')
        elif len(n_raw_h5) > 0:
            printe(' [h5r] ', red(len(n_raw_h5)), '')
        elif len(n_mtx_features) + len(n_mtx_genes) > 0:
            if len(n_mtx_features) > 0:
                printe(' [mtxz] ', red(len(n_mtx_features)), '')
            if len(n_mtx_genes) > 0:
                printe(' [mtx] ', red(len(n_mtx_genes)), '')
        elif len(n_undef_h5) > 0:
            printe(' [h5u] ', red(len(n_undef_h5)), '')

        # find the name from the catlog tree.
        def findname(root, catlog):
            for samplek in catlog['samples'].keys():
                sample = catlog['samples'][samplek]
                # print(samplek, sample['files']) # debug
                for fname in sample['files']:
                    if fname.startswith(root):
                        return {
                            'acc': sample['acc'][0],
                            'name': sample['name'][0]
                        }
                        
            snames = list(catlog['samples'].keys())
            for sname in snames:
                if root.lower().startswith(sname.lower()):
                    return {
                        'acc': catlog['samples'][sname]['acc'][0],
                        'name': catlog['samples'][sname]['name'][0]
                    }

            if root.lower().startswith(catlog['dataset']['acc'][0].lower()):
                return {
                    'acc': catlog['dataset']['acc'][0],
                    'name': 'root'
                }
            
            for fname in catlog['dataset']['files']:
                # print('dataset', catlog['dataset']['files']) # debug
                if fname.startswith(root):
                    return {
                        'acc': catlog['dataset']['acc'][0],
                        'name': 'root'
                    }
            
            return {
                'acc': 'fail',
                'name': 'fail'
            }
        
        def uniform(x):
            return x.lower() \
                .replace(' (', '(').replace(') ', ')').replace(', ',',') \
                .replace('-','.').replace('_','.') \
                .replace('   ', '.').replace('  ', '.').replace(' ', '.')
        
        def printinfo(prefix, catlog):
            fn = findname(prefix, catlog)

            printe(' ')
            
            if fn['acc'] == 'fail':
                printe(red('cannot find parent.'))
                return 0
            elif fn['acc'].lower().startswith('gse'):
                printe(formal_length(red('gse'), ':', yellow(fn['acc'].lower().replace('gse', '')), l = 16))
            elif fn['acc'].lower().startswith('gsm'):
                printe(formal_length(red('gsm'), ':', green(fn['acc'].lower().replace('gsm', '')), l = 16))
            else: printe(formal_length(red('?'), l = 16))
            
            printe(' ')

            # print name only
            if fn['acc'] == catlog['dataset']['acc'][0]:
                printe(formal_length('root', l = 30))
            elif fn['acc'] in catlog['samples'].keys():
                printe(formal_length(uniform(catlog['samples'][fn['acc']]['name'][0]), l = 30))

            printe(' ')
            
            if fn['acc'] == catlog['dataset']['acc'][0]:
                printe(formal_length('mixed' if len(catlog['dataset']['tax']) > 1 else catlog['dataset']['tax'][0], l = 6))
                if len(catlog['dataset']['tax']) > 1: return 0
                if catlog['dataset']['tax'][0] != '10090': return 0
                    
            elif fn['acc'] in catlog['samples'].keys():
                printe(formal_length(uniform(
                    'mixed' \
                        if len(catlog['samples'][fn['acc']]['tax']) > 1 \
                        else catlog['samples'][fn['acc']]['tax'][0]), l = 6))
                if len(catlog['samples'][fn['acc']]['tax']) > 1: return 0
                if catlog['samples'][fn['acc']]['tax'][0] != '10090': return 0
                
            printe(' ')

            # here, we will produce the sample summary file if necessary
            # or extract information from the chat result.
            
            atype = ''
            acc = ''
            if fn['acc'].lower().startswith('gse'):
                atype = 'gse'
                acc = fn['acc'].lower().replace('gse', '')
            elif fn['acc'].lower().startswith('gsm'):
                atype = 'gsm'
                acc = fn['acc'].lower().replace('gsm', '')
            else: return 0

            print(prefix, f'{atype}:{acc}', sep = '\t', file = _f)

            import os
            import textwrap
            
            # dump the database into raw text files.
            
            def write_sample(samp, fw):
                print(f'Sample {samp["acc"][0]}', file = fw)
                print(f'Name: {" ".join(samp["name"])}', file = fw)
                print(f'Last update: {" ".join(samp["date"])}', file = fw)
                
                parags = []
                for parag in samp["char"]:
                    parags += ['\n'.join(["    " + x for x in textwrap.wrap(parag, 76)])]
                if len(parags) == 0:
                    parags += ['    No sample characteristics available']
                print('Characteristics: \n\n' + "\n".join(parags) + '\n', file = fw)
                
                parags = []
                for parag in samp['prtcl.treat']:
                    parags += ['\n'.join(["    " + x for x in textwrap.wrap(parag, 76)])]
                if len(parags) == 0:
                    parags += ['    No treatment protocol available']
                print('Treatment protocol: \n\n' + "\n\n".join(parags) + '\n', file = fw)
                
                parags = []
                for parag in samp['prtcl.extract']:
                    parags += ['\n'.join(["    " + x for x in textwrap.wrap(parag, 76)])]
                if len(parags) == 0:
                    parags += ['    No library construction protocol available']
                print('Library protocol: \n\n' + "\n\n".join(parags) + '\n', file = fw)
                pass
                    
            if atype == 'gse':
                # the data is given completely in the gse dataset root.
                # so we need to infer sample grouping from all the samples.

                if os.path.exists(f'summary/gse/{acc}.sum'):
                    pass
                else:
                    with open(f'summary/gse/{acc}.sum', 'w') as fw:

                        # series info.
                        print(f'Series {catlog["dataset"]["acc"][0]}', file = fw)
                        print(f'Title: {" ".join(catlog["dataset"]["title"])}', file = fw)
                        print(f'Last update: {" ".join(catlog["dataset"]["date"])}', file = fw)

                        parags = []
                        for parag in catlog["dataset"]["summary"]:
                            parags += ['\n'.join(["    " + x for x in textwrap.wrap(parag, 76)])]
                        if len(parags) == 0:
                            parags += ['    No summary available']
                        print('Summary: \n\n' + "\n\n".join(parags) + '\n', file = fw)

                        parags = []
                        for parag in catlog["dataset"]["design"]:
                            parags += ['\n'.join(["    " + x for x in textwrap.wrap(parag, 76)])]
                        if len(parags) == 0:
                            parags += ['    No experimental design information available']
                        print('Design: \n\n' + "\n\n".join(parags) + '\n', file = fw)
                        print('-' * 80, file = fw)
                        
                        # samples info. (all samples)
                        for sampk in catlog['samples'].keys():
                            write_sample(catlog['samples'][sampk], fw)
                            print('-' * 80, file = fw)

                        print('Record finished.', file = fw)
                        print(f'This series contains {len(catlog["samples"])} samples.', file = fw)

            if atype == 'gsm':
                # we can just feed the exact sample information to ai.

                if os.path.exists(f'summary/gsm/{acc}.sum'):
                    pass
                else:
                    with open(f'summary/gsm/{acc}.sum', 'w') as fw:

                        # series info.
                        print(f'Series {catlog["dataset"]["acc"][0]}', file = fw)
                        print(f'Title: {" ".join(catlog["dataset"]["title"])}', file = fw)
                        print(f'Last update: {" ".join(catlog["dataset"]["date"])}', file = fw)

                        parags = []
                        for parag in catlog["dataset"]["summary"]:
                            parags += ['\n'.join(["    " + x for x in textwrap.wrap(parag, 76)])]
                        if len(parags) == 0:
                            parags += ['    No summary available']
                        print('Summary: \n\n' + "\n\n".join(parags) + '\n', file = fw)

                        parags = []
                        for parag in catlog["dataset"]["design"]:
                            parags += ['\n'.join(["    " + x for x in textwrap.wrap(parag, 76)])]
                        if len(parags) == 0:
                            parags += ['    No experimental design information available']
                        print('Design: \n\n' + "\n\n".join(parags) + '\n', file = fw)
                        print('-' * 80, file = fw)
                        
                        # samples info. (all samples)
                        write_sample(catlog['samples']['GSM' + acc], fw)
                        print('-' * 80, file = fw)

                        print('Record finished.', file = fw)
                        print(f'This series contains {len(catlog["samples"])} samples in total.', file = fw)
                        print(f'In this record, only one sample is given. But you are able to extract information with this only.', file = fw)
            
            return 1

        if len(n_filtered_h5) + len(n_raw_h5) + len(n_mtx_features) + len(n_mtx_genes) + len(n_undef_h5) == 0:
            printe(' no available data.')
            print()
            continue
        
        print()
        n_available_ds += 1
        
        if len(n_filtered_h5) > 0:
            for fp in n_filtered_h5:
                printe(formal_length('| ' + fp))
                n_available_sample += printinfo(fp + 'filtered_', fcatlog)
                print()

        elif len(n_raw_h5) > 0:
            for fp in n_raw_h5:
                printe(formal_length('| ' + fp))
                n_available_sample += printinfo(fp + 'raw_', fcatlog)
                print()

        elif len(n_mtx_features) + len(n_mtx_genes) > 0:
            if len(n_mtx_features) > 0:
                for fp in n_mtx_features:
                    printe(formal_length('| ' + fp))
                    n_available_sample += printinfo(fp + 'matrix.mtx', fcatlog)
                    print()

            if len(n_mtx_genes) > 0:
                for fp in n_mtx_genes:
                    printe(formal_length('| ' + fp))
                    n_available_sample += printinfo(fp + 'matrix.mtx', fcatlog)
                    print()

        elif len(n_undef_h5) > 0:
            for fp in n_undef_h5:
                printe(formal_length('| ' + fp))
                n_available_sample += printinfo(fp + '.h5', fcatlog)
                print()

        pass

    print()
    print(cyan('Summary: '))
    print(cyan('  Total datasets:'), red(str(n_total)))
    print(cyan('  Avaiable datasets for autoloading:'), red(str(n_available_ds)))
    print(cyan('  Available samples:'), red(str(n_available_sample)))