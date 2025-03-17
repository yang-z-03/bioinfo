
import os
from multiprocessing import Pool as pool
from time import sleep
from printtbl import printtbl
from ansi import common_length
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import font_manager as fm
import matplotlib as mpl
import pandas as pd
import datetime

ftpath = "/home/data/yangzhen/fonts/arial.ttf"
ftprop = fm.FontProperties(fname = ftpath)
rc = {'font.sans-serif': ['Arial', 'DejaVu Sans', 'Bitstream Vera Sans']}
sns.set(font = 'sans-serif', rc = rc)

_fe = open('configs/qc.err', 'w')
_fe.writelines([f'[start] {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")} \n\n'])
_fe.flush()

def status_qc(registry):

    # get the download status as a string array
    # which will be appended to the list display.

    stat = []
    src = os.listdir('processed/qc')
    adatas = os.listdir('adata')
    for i in range(len(registry['name'])):
        name = registry['name'][i]
        if name + '.h5ad' in src:
            stat += ['done']
            continue
        elif name + '.h5ad' in adatas:
            stat += ['pending']
        else: stat += ['.']
    
    return stat

def styles_qc(x: str):
    if x == 'done':
        return '\033[32;1m' + x + '\033[0m'
    elif x == 'pending':
        return '\033[33;1m' + x + '\033[0m'
    else: return '\033[31;0m' + x + '\033[0m'


def display(registry, n_row, show_status = False):

    from metagen import status_dtype, styles_dtype
    from download import status_data, styles_data

    registry['dtype'] = status_dtype(registry)
    registry['data'] = status_data(registry)
    registry['qc'] = status_qc(registry)
    from list import styles_acctype, styles_number
    total_row = len(registry['data'])
    nfrom = 0

    while nfrom < total_row:

        cols = printtbl(
            registry, styles = { 
                'type': styles_acctype, 'accession': styles_number,
                'dtype': styles_dtype, 'data': styles_data, 'qc': styles_qc
            },
            appends = { 
                'type': 9, 'accession': 11,
                'dtype': 11, 'data': 11, 'qc': 11
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
        def gprintc(x, ln, l = 40):
            ansi_move_cursor(ln - 1, 0)
            rprintc(x, l)
            ansi_move_cursor(1 - ln, 0)
        
        qcs = registry['qc'][nfrom : min(total_row, nfrom + n_row)]
        running_task = []
        running_id = []
        qclogs = []
        for s_name, s_acctype, s_acc, s_data, s_dtype, s_qc, s_n in \
            zip(
                registry['name'][nfrom : min(total_row, nfrom + n_row)], 
                registry['type'][nfrom : min(total_row, nfrom + n_row)], 
                registry['accession'][nfrom : min(total_row, nfrom + n_row)], 
                registry['data'][nfrom : min(total_row, nfrom + n_row)], 
                registry['dtype'][nfrom : min(total_row, nfrom + n_row)],
                registry['qc'][nfrom : min(total_row, nfrom + n_row)],
                range(len(qcs))
            ):

            if s_qc == 'done':
                gprintc('done.', s_n + 1)
                continue
            elif s_qc == 'pending':
                gprintc(f'running {s_name} ...', s_n + 1)
                running_task += [f'{s_name}.h5ad']
                running_id += [s_n + 1]
                qclogs += [f'processed/logs/qc/{s_name}.qclog']
                continue
            else: gprintc('skipped.', s_n + 1)
            pass

        # set up the thread pool for simutaneous qc.
        with pool(25) as p:

            # using thread pool
            results = [p.apply_async(qc, args = (fn, id, qcl)) 
                       for fn, id, qcl in zip(running_task, running_id, qclogs)]
            final = [res.get() for res in results]

            # do not. for debug only
            # final = []
            # for fn, id, qcl in zip(running_task, running_id, qclogs):
            #     final += [qc(fn, id, qcl)]
            
            for f in final: gprintc('done.', f)
        
        ansi_move_cursor(n_row, 0)
        print('\r', end = '\n')

        nfrom += n_row

    pass

def qc(fname, run_id, qclog, run_scrublet = True):

    from contextlib import redirect_stdout

    with open(qclog, 'w') as qclogf:
        with redirect_stdout(qclogf):

            # i cannot find a way to shut up omicverse ...
            import omicverse as ov
            import scanpy as sc
            import sys
            import warnings
            import numpy as np
            import matplotlib.pyplot as plt
            import anndata as ad
            import pickle
            from matplotlib.ticker import FuncFormatter as fformat
            warnings.filterwarnings('ignore')
            sc.settings.verbosity = 0

            name = fname.replace('.h5ad', '')

            adata = sc.read('adata/' + fname)

            print('{0} samples, {1} genes. total {2}'.format(
                adata.n_obs, adata.n_vars, adata.n_obs * adata.n_vars
            ))

            print('{0} nan\'s introduced by coercing. ({0}/{1})'.format(
                np.isnan(adata.X.data).sum(), adata.n_obs * adata.n_vars
            ))
            
            # if the anndata goes with too small gene detection (commonly due to name mismatch
            # and resulting in zero genes,) we will just skip this dataset.
            if adata.n_vars <= 200 or adata.n_obs <= 400:
                print('skipping the dataset for quality control due to insufficient cell numbers or genes.')
                _fe.writelines([f'[{name.lower()}] skipping the dataset for quality control due to insufficient cell numbers or genes.' + '\n'])
                _fe.flush()
                return run_id

            # the coercing of sample merging introduced nan's in the data frame.
            # we need to replace them with 0 before performing later steps.
            adata.X.data = np.nan_to_num(adata.X.data, copy = False)

            # we should extract ribosomal genes and mitochondrial genes from g***.

            mito_genes = [
                "g17895", "g17896", "g17897", "g17898", "g17899", "g17900", 
                "g17901", "g17902", "g17903", "g17904", "g17905", "g17906", 
                "g17907", "g17908", "g17909", "g17910", "g17911", "g17912", 
                "g17913", "g17914", "g17915", "g17916", "g17917", "g17918", 
                "g17919", "g17920", "g17921", "g17922", "g17923", "g17924", 
                "g17925", "g17926", "g17927", "g17928", "g17929", "g17930", 
                "g18691"
            ]

            ribo_genes = [
                "g127", "g439", "g583", "g617", "g652", "g979", "g1137", "g1183",
                "g1184", "g1185", "g1280", "g1395", "g1409", "g1798", "g2067", 
                "g2531", "g2560", "g2952", "g4153", "g4355", "g4581", "g4699", 
                "g4969", "g4996", "g5082", "g5144", "g5193", "g5250", "g6477",
                "g6562", "g6714", "g6980", "g7402", "g7826", "g7991", "g8271", 
                "g8276", "g8936", "g9009", "g9841", "g10302", "g10446", "g10753", 
                "g10842", "g11319", "g11352", "g11447", "g11454", "g12099", 
                "g12240", "g12374", "g12925", "g13162", "g13312", "g13336", 
                "g13662", "g13671", "g13683", "g13739", "g13770", "g13914", 
                "g14043", "g14044", "g14125", "g14340", "g14528", "g14602", 
                "g14824", "g15079", "g15122", "g15210", "g15435", "g15857", 
                "g16213", "g16214", "g16234", "g16328", "g16337", "g16406", 
                "g16413", "g16446", "g16450", "g16566", "g16570", "g16590", 
                "g16592", "g16604", "g16605", "g16606", "g16629", "g16661",
                "g16700", "g16709", "g16726", "g16759", "g16865", "g16870", 
                "g16871", "g16884", "g16911", "g16936", "g17007", "g17017", 
                "g17021", "g17035", "g17051", "g17065", "g17130", "g17173", 
                "g17230", "g17234", "g17274", "g17275", "g17295", "g17339", 
                "g17341", "g17344", "g17369", "g17386", "g17432", "g17457", 
                "g17459", "g17491", "g17501", "g17509", "g17559", "g17582", 
                "g17616", "g17652", "g17685", "g17781", "g17787", "g17878", 
                "g18692", "g18774", "g18782", "g18799", "g18809", "g18825", 
                "g18885", "g18901", "g18904", "g18914", "g18941", "g18945", 
                "g18953", "g18954", "g18988", "g19043", "g19155", "g19180", 
                "g19197", "g19258", "g19279", "g19281", "g19405", "g19425", 
                "g19616", "g19623", "g19643", "g19687", "g19723", "g19740", 
                "g19783", "g19806", "g19807", "g19817", "g19944", "g20011", 
                "g20017", "g20261", "g20277", "g20449", "g20458", "g20614", 
                "g20976", "g21862", "g21870", "g22260", "g22369", "g22471", 
                "g22521", "g22530", "g22535", "g22538", "g22555", "g22677", 
                "g22689", "g22697", "g22786", "g22788", "g22800", "g22804", 
                "g22808", "g22838", "g22948", "g22973", "g22984", "g22989", 
                "g22997", "g23016", "g23095", "g23157", "g23191", "g23195", 
                "g23241", "g23244", "g23263", "g23281", "g23296", "g23309", 
                "g23373", "g23396", "g23446", "g23454", "g23467", "g23507", 
                "g23518", "g23535", "g23564", "g23653", "g23692", "g23701", 
                "g23718", "g23721", "g23757", "g23767", "g23815", "g23841", 
                "g23880", "g23913", "g23935", "g23955", "g24190", "g24201", 
                "g24226", "g24245", "g24328", "g24333", "g24345", "g24410", 
                "g24424", "g24437", "g24438", "g24465", "g24466", "g24486", 
                "g24487", "g24548", "g24565", "g24571", "g24623", "g24636", 
                "g24643", "g24647", "g24673", "g24752", "g24937", "g25015", 
                "g25027", "g25041", "g25076", "g25158", "g25170", "g25282", 
                "g25309", "g25333", "g25337", "g25447", "g25515", "g25539", 
                "g25637", "g25687", "g25688", "g25753", "g25802", "g25830", 
                "g25834", "g25855", "g25872", "g25881", "g25891", "g26773", 
                "g26833", "g29292", "g29326", "g29343", "g29442", "g29501", 
                "g29599", "g29726", "g29738", "g29829", "g29840", "g29882", 
                "g29969", "g30270", "g30288", "g30297", "g30344", "g30633", 
                "g30668", "g30981", "g31142", "g31219", "g31260", "g31285", 
                "g31324", "g31335", "g31362", "g31389", "g31424", "g31451", 
                "g31461", "g31605", "g31660", "g31720", "g31823", "g31829", 
                "g31842", "g31896", "g31902", "g31964", "g32039", "g32074", 
                "g32129", "g32403", "g32463", "g32467", "g32509", "g32546", 
                "g32669", "g32681", "g32724", "g32758", "g32811", "g32824", 
                "g32930", "g33023", "g33029", "g33057", "g33084", "g33218", 
                "g33304", "g33333", "g33386", "g33444", "g33703", "g33783", 
                "g34040", "g34153", "g34780", "g35488", "g35697", "g35716", 
                "g36080", "g36260", "g36261", "g36425", "g36426", "g36932", 
                "g37002", "g37021", "g37096", "g37455", "g38640", "g38679", 
                "g38928", "g38991", "g39270", "g39370", "g39675", "g39760", 
                "g39952", "g39976", "g40038", "g40751", "g40772", "g40949", 
                "g41317", "g42205", "g42641", "g42679", "g42725", "g43020", 
                "g43690", "g44596", "g44837", "g45263", "g45369", "g45416", 
                "g45487", "g45841", "g46080", "g46553", "g46739", "g47031", 
                "g47045", "g47265", "g47903", "g47914", "g47918", "g48199", 
                "g48283", "g48647", "g49450", "g49532", "g49567", "g49672", 
                "g49902", "g49963", "g50218", "g50236", "g50409", "g50776", 
                "g50972", "g51005", "g51030", "g51111", "g51145", "g51276", 
                "g51317", "g51718", "g51821", "g52104", "g52214"
            ]

            # we should remove a sample if it contains less than 200 cells.
            # this should not be normal for a typical 10x like platform.
            cell_counts = adata.obs['sample'].value_counts()
            sample_names = cell_counts.index.tolist()
            sample_size = cell_counts.values.tolist()
            smaller = []
            
            for nm, sz in zip(sample_names, sample_size):
                if sz < 200: smaller.append(nm)

            if len(smaller) > 0:
                sname = adata.obs['sample'].tolist()
                keep = [(x not in smaller) for x in sname]
                adata = adata[keep, :].copy() # remove samples with too small size.
                print('remove sample', smaller, 'due to <200 cells')
                print('remove sample', smaller, 'due to <200 cells')
                _fe.writelines([f'[{name.lower()}] remove sample {smaller} due to <200 cells.' + '\n'])
                _fe.flush()
                
            adata.var['is.mito'] = [x in mito_genes for x in adata.var_names.tolist()]
            adata.var['is.ribo'] = [x in ribo_genes for x in adata.var_names.tolist()]

            print('found mitochondria genes: {}'.format(adata.var['is.mito'].sum()))
            print('found ribosomal genes: {}'.format(adata.var['is.ribo'].sum()))

            sample_metrics, gene_metrics = sc.pp.calculate_qc_metrics(
                adata, qc_vars = ['is.mito', 'is.ribo'],
                inplace = False, percent_top = [20], log1p = False
            )

            adata.obs['pct.mito'] = sample_metrics['pct_counts_is.mito']
            adata.obs['pct.ribo'] = sample_metrics['pct_counts_is.ribo']
            adata.obs['n.umi'] = sample_metrics['total_counts']
            adata.obs['n.det'] = sample_metrics['n_genes_by_counts']
            adata.obs['pct.top.20'] = sample_metrics['pct_counts_in_top_20_genes']

            if run_scrublet:
                try:
                    sc.external.pp.scrublet(adata, random_state = 42, batch_key = 'sample')
                    adata.obs.rename(
                        columns = {'doublet_score': 'doublet.score', 'predicted_doublet': 'is.doublet'}, 
                        inplace = True
                    )
                
                except:
                    _fe.writelines([f'[{name.lower()}] error running scrublet. \n'])
                    _fe.flush()
                    adata.obs['doublet.score'] = 0
                    adata.obs['is.doublet'] = False
            else:
                adata.obs['doublet.score'] = 0
                adata.obs['is.doublet'] = False

            adata.obs['pass.qc.mito'] = False
            adata.obs['pass.qc.umi'] = False
            adata.obs['pass.qc.detection'] = False

            samples = adata.obs['sample'].value_counts().index.tolist()

            # split the whole dataset into samples. and calculate the hvgs separately,
            # as well as to define quality control criteria separately. each sample
            # in the dataset will have a pair of quality map plotted.

            thresholds = {
                'mt.up': 0.20 * 100,
                'n.umi': 400,
                'n.det': 200
            }

            def format_tick_labels(x, pos):
               return '{0:.1f}k'.format(x / 1000)

            # change the color set to category.
            adata.obs['is.doublet'] = adata.obs['is.doublet'].astype('category')
            splices = []
            var_table = adata.var
            
            for s in samples:

                print('')
                splice = adata[adata.obs['sample'] == s,:]
                print('splicing [{0}]: {1} samples, {2} genes. total {3}'.format(
                    s, splice.n_obs, splice.n_vars, splice.n_obs * splice.n_vars
                ))

                splice.obs['pass.qc.mito'] = splice.obs['pct.mito'] < thresholds['mt.up']
                splice.obs['pass.qc.umi'] = ov.pp._qc.mads_test(splice.obs, 'n.umi', nmads = 6, lt = thresholds)
                splice.obs['pass.qc.detection'] = ov.pp._qc.mads_test(splice.obs, 'n.det', nmads = 6, lt = thresholds)  

                counts_thr = ov.pp._qc.mads(splice.obs, 'n.umi', nmads = 6, lt = thresholds)
                dets_thr = ov.pp._qc.mads(splice.obs, 'n.det', nmads = 6, lt = thresholds)

                n1 = splice.shape[0]

                print('[{}] criteria: counts: ({}, {}); keeping {} / {}'.format(
                    s, counts_thr[0], counts_thr[1], np.sum(splice.obs['pass.qc.umi']), n1
                ))

                print('[{}] criteria: genes: ({}, {}); keeping {} / {}'.format(
                    s, dets_thr[0], dets_thr[1], np.sum(splice.obs['pass.qc.detection']), n1
                ))

                print('[{}] criteria: mito percentage: {}; keeping {} / {}'.format(
                    s, thresholds['mt.up'], np.sum(splice.obs['pass.qc.mito']), n1
                ))

                splice.obs['qc'] = (splice.obs['pass.qc.mito']) & (splice.obs['pass.qc.umi']) & (splice.obs['pass.qc.detection'])

                splices += [splice]

                # plot the essential figures for visualization:

                dataframe = {
                    'n.umi': splice.obs['n.umi'].tolist(),
                    'n.det': splice.obs['n.det'].tolist(),
                    'pct.mito': splice.obs['pct.mito'].tolist(),
                    'is.doublet': splice.obs['is.doublet'].tolist()
                }

                dataframe = pd.DataFrame(dataframe)
                dataframe['is.doublet'] = dataframe['is.doublet'].astype('category')

                if not os.path.exists('processed/graphs/mito/{0}'.format(name.lower())):
                    os.makedirs('processed/graphs/mito/{0}'.format(name.lower()))

                if not os.path.exists('processed/graphs/detections/{0}'.format(name.lower())):
                    os.makedirs('processed/graphs/detections/{0}'.format(name.lower()))

                try:
                    plot1 = sns.jointplot(data = dataframe, x = 'n.umi', y = 'pct.mito', hue = 'is.doublet', 
                                          legend = False, palette = ['black', 'red'], height = 4, s = 10)
                    # plot1.plot_joint(sns.kdeplot, color = 'r', zorder = 99, levels = 6)
                    plot1.refline(y = thresholds['mt.up'], color = 'crimson')
                    plot1.refline(x = counts_thr[0], color = 'slateblue')
                    plot1.refline(x = counts_thr[1], color = 'rebeccapurple')
                    plt.xscale('log')
                    plt.yscale('log')
                    plt.xticks(fontproperties = ftprop); plt.yticks(fontproperties = ftprop)
                    plt.title(''); 
                    plt.xlabel('UMI detection', fontproperties = ftprop)
                    plt.ylabel('Mitochondrial ratio', fontproperties = ftprop)
                    plt.savefig('processed/graphs/mito/{0}/{1}.png'.format(name.lower(), s.replace(':', '-')))
                except:
                    _fe.writelines([f'[{name.lower()}/{s}] error when creating mito.png. \n'])
                    _fe.flush()
                    pass
                    
                try:
                    plot2 = sns.jointplot(data = dataframe, x = 'n.umi', y = 'n.det', hue = 'is.doublet',
                                          legend = False, palette = ['black', 'red'], height = 4, s = 10)
                    # plot2.plot_joint(sns.kdeplot, color = 'r', zorder = 99, levels = 6)
                    plot2.refline(y = dets_thr[0], color = 'crimson')
                    plot2.refline(y = dets_thr[1], color = 'orangered')
                    plot2.refline(x = counts_thr[0], color = 'slateblue')
                    plot2.refline(x = counts_thr[1], color = 'rebeccapurple')
                    plt.xscale('log')
                    plt.yscale('log')
                    plt.xticks(fontproperties = ftprop); plt.yticks(fontproperties = ftprop)
                    plt.title('')
                    plt.xlabel('UMI detection', fontproperties = ftprop)
                    plt.ylabel('Gene detection', fontproperties = ftprop)
                    plt.savefig('processed/graphs/detections/{0}/{1}.png'.format(name.lower(), s.replace(':', '-')))
                except:
                    _fe.writelines([f'[{name.lower()}/{s}] error when creating detections.png. \n'])
                    _fe.flush()
                    pass
                    
                pass

            merged = ad.concat(splices)
            # here, we ensure that the splicing operation did not change the order
            # of genes. so we can paste the gene table per se.
            # while the sample ordering is shuffled.
            merged.var = var_table

            if run_scrublet:
                if 'scrublet' in adata.uns.keys(): del adata.uns['scrublet']
                if 'n_genes' in adata.obs.keys(): del adata.obs['n_genes']

            # by now, before any modification on subset operation, we will save the
            # object. combine these criteria into one.

            merged.obs.to_csv('processed/metrics/{0}.tsv'.format(name), sep = '\t', index = False)
            qc = merged[merged.obs['qc'],:]
            qc.write_h5ad('processed/qc/{0}.h5ad'.format(name))

            # calculate hvg and make the union.
            # this is sometimes optional.
            hvgs = []
            for i in range(len(samples)):

                print('')

                s = samples[i]
                spl = splices[i]

                n2, g0 = spl.shape
                sc.pp.filter_cells(spl, min_genes = 200)
                n3, g0 = spl.shape
                sc.pp.filter_genes(spl, min_cells = 3)
                n4, g1 = spl.shape

                print(f'[{s}] gene filter: (from) c{n2}:g{g0} (keeping) c{n4}:g{g1}')

                spl.layers['counts'] = spl.X

                try:
                    seurat_v3_hvg = sc.pp.highly_variable_genes(
                        spl,
                        flavor = 'seurat_v3',
                        layer = 'counts',
                        n_top_genes = 2000,
                        subset = False,
                        inplace = False,
                    )

                    varnames = spl.var_names.copy()
                    ls = set(varnames[seurat_v3_hvg['highly_variable'] == True].tolist())
                    hvgs += [ls]
                    print(f'[{s}] {len(ls)} hvgs found (by seurat vst).')

                except: print(f'[{s}] hvg picking failed.')

            hvgd = []
            for l in hvgs: hvgd += list(l)
            hvg = set(hvgd)

            with open('processed/hvg/{0}.pkl'.format(name), 'wb') as f:
                pickle.dump(hvg, f)

            return run_id
