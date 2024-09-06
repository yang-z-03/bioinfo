
def scrublet_py(i, j, val, dim, expected_doublet_rate, min_counts, 
				min_cells, min_gene_variability_pctl, n_prin_comps,
				sim_doublet_ratio, n_neighbors):
	
	# print('<type> i : ', type(i))
	# print('<type> j : ', type(j))
	# print('<type> val : ', type(val))
	# print('<type> dim : ', type(dim))
	# print('<type> expected_doublet_rate : ', type(expected_doublet_rate))
	# print('<type> min_counts : ', type(min_counts))
	# print('<type> min_cells : ', type(min_cells))
	# print('<type> min_gene_variability_pctl : ', type(min_gene_variability_pctl))
	# print('<type> n_prin_comps : ', type(n_prin_comps))
	# print('<type> sim_doublet_ratio : ', type(sim_doublet_ratio))
	# print('<type> n_neighbors : ', type(n_neighbors))

	# import pickle

	# with open('scrublet-dump', 'wb') as f:
	# 	pickle.dump(
	# 		[i, j, val, dim, expected_doublet_rate, min_counts, 
	# 		 min_cells, min_gene_variability_pctl, n_prin_comps,
	# 		 sim_doublet_ratio, n_neighbors], f
	# 	)
	
	import matplotlib
	matplotlib.use('agg')
	import scrublet as scr
	import scipy.io
	import numpy as np
	import os
	from scipy.sparse import csc_matrix
	
	data = csc_matrix((val, (i, j)), shape = dim)
	scrub = scr.Scrublet(data, expected_doublet_rate = expected_doublet_rate,
					     sim_doublet_ratio = int(sim_doublet_ratio),
						 n_neighbors = int(n_neighbors))
	
	doublet_scores, predicted_doublets = scrub.scrub_doublets(
		min_counts = int(min_counts),
		min_cells = int(min_cells),
		min_gene_variability_pctl = min_gene_variability_pctl,
		n_prin_comps = int(n_prin_comps)
	)
	
	return (doublet_scores, predicted_doublets)
