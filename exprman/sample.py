
# produce a sample table with linked files.
# we will not use the printtbl routines, but rather print by ourself.

import os
import pandas as pd

_table_prev = None
if os.path.exists('configs/lookup.tsv'):
    _table_prev = pd.read_table('configs/lookup.tsv')
    _table_prev['prefix'] = _table_prev['prefix'].astype('str')
    os.remove('configs/lookup.tsv')
    
_f = open('configs/lookup.tsv', 'w')

print(
    'accession\tdataset\tsample\tdstype\ttax\t'
    'dtype\tfiltered\tselected\tencode\tname\tprefix\t'
    'strain\tsex\ttissue\tage\texpm\tgenetic\tsort',
    file = _f
)

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

def load_sample_config():
    
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

    return config


def display(registry, n_row, show_status = False):

    n_total = 0
    n_available_sample = 0
    manual_config = load_sample_config()
    
    for s_name, s_type, s_acc in zip(
        registry['name'], 
        registry['type'], 
        registry['accession']
    ):

        if s_type != 'gse':
            # print nothing for unsupported dataset. just skip it.
            # printe(formal_length(red(s_type), ':', yellow(s_acc), l = 40))
            # printe(formal_length(' dataset not supported.'))
            # print()
            continue

        gserec = 'GSE' + s_acc
        import GEOparse as geo
        gse = geo.get_GEO(gserec, destdir = 'configs/geo', silent = True)

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

        if not os.path.exists(f'src/{s_name}'):
            continue

        n_total += 1
        # list downloaded files
        down = os.listdir(f'src/{s_name}')
        visible = []
        for x in down:
            if not x.startswith('.'): visible += [x]
        
        n_mtx_genes = []
        n_mtx_features = []
        n_h5 = []
        unrecog = []
        removed = []

        # filter and test the name conventions.
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

        # precedency rules
        # h5f > h5r > mtx + mtxz （filtered > raw) > h5u

        # if sample contains its own set of matrix, we will ignore those provided
        # as a whole in the dataset root. if all samples provide no valid matrix
        # we will search the root directory.

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
        
        def getinfo(prefix, catlog, manual):
            fn = findname(prefix, catlog)
            
            metadata = {
                'accession': '',  # gse series
                'dataset': s_name.replace(str(s_acc) + '-', '').replace(str(s_acc), ''), 
                'sample': '',     # gsm sample, or 'root' if unknown.
                'name': '',       # original sample name given by the database
                'tax': '',        # taxonomy code, or mixed if more than one.
                'dtype': '',      # file format, mtx, mtxz or h5
                'dstype': '',     # dataset type, scrna etc.
                # 'filtered': '', # if this is filtered. raw, filt or ?
                # 'selected': '', # if this should be selected automatically. x or .
                                  # note that all files in available extensions will be listed in this table
                                  # but only some of the de-duplicated ones may be processed next.
                                  # this is selected using our precedence rule.
                
                'prefix': '',     # the file prefix for selection (without the trailing matrix.tsv[.gz] or .h5[.gz]
                # 'encode': '',   # unique sample encoding
                'strain': '?',    # strain
                'sex': '?',       # sex code, m for male, f for female, ? for unknown
                'tissue': '?',    # tissue encode
                'age': '?',       # age
                'expm': '?',      # experiment grouping
                'genetic': '?',
                'sort': '?'
            }
            
            if fn['acc'] == 'fail':
                
                # this file exists, but with no explicit evidence (without inferrence) that
                # it belongs to any of the samples.  in this case, we should only provide
                # metadata information by ourself in our hand-make samples file.
                
                metadata['accession'] = 'gse:' + catlog['dataset']['acc'][0][3:]
                metadata['sample'] = 'root'
                fn['acc'] = catlog['dataset']['acc'][0]
                fn['name'] = 'root'
                
            elif fn['acc'].lower().startswith('gse'):
                metadata['accession'] = 'gse:' + catlog['dataset']['acc'][0][3:]
                metadata['sample'] = 'root'
                
            elif fn['acc'].lower().startswith('gsm'):
                metadata['accession'] = 'gse:' + catlog['dataset']['acc'][0][3:]
                metadata['sample'] = 'gsm:' + fn['acc'][4:]
            
            else:
                metadata['accession'] = 'gse:' + catlog['dataset']['acc'][0][3:]
                metadata['sample'] = 'root'
            
            # print name only
            # you should only set the name manually!
                
            # if fn['acc'] == catlog['dataset']['acc'][0]:
            #     metadata['name'] = 'root'
            # elif fn['acc'] in catlog['samples'].keys():
            #     metadata['name'] = catlog['samples'][fn['acc']]['name'][0]
            # else: metadata['name'] = 'root'

            if fn['acc'] == catlog['dataset']['acc'][0]:
                if len(catlog['dataset']['tax']) > 1: metadata['tax'] = 'mixed'
                else: metadata['tax'] = catlog['dataset']['tax'][0]
                    
            elif fn['acc'] in catlog['samples'].keys():
                if len(catlog['samples'][fn['acc']]['tax']) > 1: metadata['tax'] = 'mixed'
                else: metadata['tax'] = catlog['samples'][fn['acc']]['tax'][0]

            # extract sample information.
            metadata['prefix'] = prefix
            
            # the manual configuration overwrites the machine-generated ones.
            def get_props(s_name, s_sample, config):
                
                props_items = ['name', 'tissue', 'genetic', 'strain', 'sort', 'expm', 'age', 'sex']
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
                
                pdict['dataset'] = s_name.replace(str(s_acc) + '-', '').replace(str(s_acc), '')
                return pdict, error_code
            
            # the machine-generated configs is stored in query/gse(gsm)/***.json, which
            # is the json-formatted metadata given by llm ai model.
                
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
            
            else: # only check if manually given, then return
                prop, e = get_props(s_name, prefix, manual)
                if e != 0: return metadata
                if 'strain' in prop.keys(): metadata['strain'] = prop['strain']
                if 'sex' in prop.keys(): metadata['sex'] = prop['sex']
                if 'tissue' in prop.keys(): metadata['tissue'] = prop['tissue']
                if 'age' in prop.keys(): metadata['age'] = prop['age']
                if 'expm' in prop.keys(): metadata['expm'] = prop['expm']
                if 'sort' in prop.keys(): metadata['sort'] = prop['sort']
                if 'genetic' in prop.keys(): metadata['genetic'] = prop['genetic']
                if 'name' in prop.keys(): metadata['name'] = prop['name']
                return metadata

            import json

            # set dstype
            if os.path.exists(f'configs/dtype/{metadata["accession"].replace(":", "-")}.json'):
                with open(f'configs/dtype/{metadata["accession"].replace(":", "-")}.json', 'r', encoding = 'utf-8') as fj:
                    try:
                        jobject = json.loads(fj.read())['samples']

                        if atype == 'gsm':
                            for xk in jobject.keys():
                                if f'gsm{acc}' == xk.lower():
                                    metadata['dstype'] = jobject[xk]
                        elif atype == 'gse':
                            evalstat = []
                            for xk in jobject.keys():
                                if jobject[xk] not in evalstat: evalstat += [jobject[xk]]
                            metadata['dstype'] = ' '.join(evalstat).lower()
                            
                    except: pass

            # set tissue
            if os.path.exists(f'configs/tissue/{metadata["accession"].replace(":", "-")}.json'):
                with open(f'configs/tissue/{metadata["accession"].replace(":", "-")}.json', 'r', encoding = 'utf-8') as fj:
                    try:
                        jobject = json.loads(fj.read())['samples']

                        if atype == 'gsm':
                            for xk in jobject.keys():
                                if f'gsm{acc}' == xk.lower():
                                    metadata['tissue'] = jobject[xk]
                        elif atype == 'gse':
                            evalstat = []
                            for xk in jobject.keys():
                                if jobject[xk].replace(' ', '.') not in evalstat: evalstat += [jobject[xk].replace(' ', '.')]
                            metadata['tissue'] = ' '.join(evalstat).lower()
                            
                    except: pass

            # set strain
            if os.path.exists(f'configs/strain/{metadata["accession"].replace(":", "-")}.json'):
                with open(f'configs/strain/{metadata["accession"].replace(":", "-")}.json', 'r', encoding = 'utf-8') as fj:
                    try:
                        jobject = json.loads(fj.read())['samples']

                        if atype == 'gsm':
                            for xk in jobject.keys():
                                if f'gsm{acc}' == xk.lower():
                                    metadata['strain'] = jobject[xk]
                        elif atype == 'gse':
                            evalstat = []
                            for xk in jobject.keys():
                                if jobject[xk].replace(' ', '.') not in evalstat: evalstat += [jobject[xk].replace(' ', '.')]
                            metadata['strain'] = ' '.join(evalstat).lower()
                            
                    except: pass
                        
            prop, e = get_props(s_name, prefix, manual)
            if e == 0:
                if 'strain' in prop.keys(): metadata['strain'] = prop['strain']
                if 'sex' in prop.keys(): metadata['sex'] = prop['sex']
                if 'tissue' in prop.keys(): metadata['tissue'] = prop['tissue']
                if 'age' in prop.keys(): metadata['age'] = prop['age']
                if 'expm' in prop.keys(): metadata['expm'] = prop['expm']
                if 'sort' in prop.keys(): metadata['sort'] = prop['sort']
                if 'genetic' in prop.keys(): metadata['genetic'] = prop['genetic']
                if 'name' in prop.keys(): metadata['name'] = prop['name']

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
            
            return metadata

        def printinfo(metas):
            if len(metas) == 0: return
            _acc = _table_prev.loc[_table_prev['accession'] == metas[0]['accession'], :]
            sample_id = len(_acc)
            for m in metas:
                
                m2 = m.copy()
                if _table_prev is not None:
                    _row = _acc.loc[_acc['prefix'] == (m['prefix'] if len(m['prefix']) > 0 else 'nan'), :]
                    if len(_row) == 1:
                        m2['encode'] = _row['encode'].tolist()[0]
                    
                    else: 
                        sample_id += 1
                        m2['encode'] = 's:' + str(sample_id)
                        print(f'find new sample {m2["encode"]} for [{m["accession"]}/{m["prefix"]}')
                
                else: 
                    sample_id += 1
                    m2['encode'] = 's:' + str(sample_id)
                    print(f'find new sample {m2["encode"]} for [{m["accession"]}/{m["prefix"]}')
                
                fmt = '{accession}\t{dataset}\t{sample}\t{dstype}\t{tax}\t' \
                      '{dtype}\t{filtered}\t{selected}\t{encode}\t{name}\t{prefix}\t' \
                      '{strain}\t{sex}\t{tissue}\t{age}\t{expm}\t{genetic}\t{sort}'
                
                print(fmt.format(**m2), file = _f)

        # not even one available data for autoloading.
        if len(n_h5) + len(n_mtx_features) + len(n_mtx_genes) == 0:
            continue

        metas = []
        id_root = 1
        
        h5r = 0; h5f = 0; h5o = 0
        mtxzr = 0; mtxzf = 0; mtxzo = 0
        mtxr = 0; mtxf = 0; mtxo = 0
        
        if len(n_h5) > 0:
            for fp in n_h5:
                meta = getinfo(fp, fcatlog, manual_config)
                meta['encode'] = f's:{id_root}'
                id_root += 1
                meta['dtype'] = 'h5'
                if 'raw_' in fp:
                    meta['filtered'] = 'raw'
                    h5r += 1
                elif 'filtered_' in fp:
                    meta['filtered'] = 'filt'
                    h5f += 1
                else: 
                    meta['filtered'] = '?'
                    h5o += 1
                metas += [meta]

        if len(n_mtx_features) > 0:
            for fp in n_mtx_features:
                meta = getinfo(fp, fcatlog, manual_config)
                meta['encode'] = f's:{id_root}'
                id_root += 1
                meta['dtype'] = 'mtxz'
                if 'raw_' in fp: 
                    meta['filtered'] = 'raw'
                    mtxzr += 1
                elif 'filtered_' in fp:
                    meta['filtered'] = 'filt'
                    mtxzf += 1
                else:
                    meta['filtered'] = '?'
                    mtxzo += 1
                metas += [meta]

        if len(n_mtx_genes) > 0:
            for fp in n_mtx_genes:
                meta = getinfo(fp, fcatlog, manual_config)
                meta['encode'] = f's:{id_root}'
                id_root += 1
                meta['dtype'] = 'mtx'
                if 'raw_' in fp:
                    meta['filtered'] = 'raw'
                    mtxr += 1
                elif 'filtered_' in fp:
                    meta['filtered'] = 'filt'
                    mtxf += 1
                else: 
                    meta['filtered'] = '?'
                    mtxo += 1
                metas += [meta]

        # finally, check the list of sample, and apply the precedence rule on which
        # should be taken into the atlas.
        # precedency rules
        # h5f > h5r > mtx + mtxz （filtered > raw) > h5u

        for i in range(len(metas)):
            metas[i]['selected'] = '.'
        
        if h5f > 0:
            for i in range(len(metas)):
                if metas[i]['dtype'] == 'h5' and metas[i]['filtered'] == 'filt':
                    metas[i]['selected'] = 'x'
        
        elif h5r > 0:
            for i in range(len(metas)):
                if metas[i]['dtype'] == 'h5' and metas[i]['filtered'] == 'raw':
                    metas[i]['selected'] = 'x'
                    
        elif mtxzf + mtxf > 0:
            for i in range(len(metas)):
                if metas[i]['dtype'] in ['mtxz', 'mtx'] and metas[i]['filtered'] == 'filt':
                    metas[i]['selected'] = 'x'
                    
        elif mtxzr + mtxr > 0:
            for i in range(len(metas)):
                if metas[i]['dtype'] in ['mtxz', 'mtx'] and metas[i]['filtered'] == 'raw':
                    metas[i]['selected'] = 'x'

        elif h5o > 0:
            for i in range(len(metas)):
                if metas[i]['dtype'] == 'h5' and metas[i]['filtered'] == '?':
                    metas[i]['selected'] = 'x'

        elif mtxzo + mtxo > 0:
            for i in range(len(metas)):
                if metas[i]['dtype'] in ['mtxz', 'mtx'] and metas[i]['filtered'] == '?':
                    metas[i]['selected'] = 'x'

        # only murine dataset.
        for i in range(len(metas)):
            if metas[i]['tax'] != '10090':
                metas[i]['selected'] = '.'

        printinfo(metas)
        for i in range(len(metas)):
            if metas[i]['selected'] == 'x':
                n_available_sample += 1
        pass

    print(cyan('Summary: '))
    print(cyan('  Total datasets:'), red(str(n_total)))
    print(cyan('  Available samples:'), red(str(n_available_sample)))
    