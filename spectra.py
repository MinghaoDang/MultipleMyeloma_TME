import Spectra
import pandas as pd
import scanpy as sc
from collections import Counter
from Spectra import K_est as kst
import pickle
import re
#
#adata = sc.AnnData(X = count)
adata = sc.read_h5ad('MM.TME.235samples.TME.no_Eryth.h5ad')
adata.var_names_make_unique()
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=30)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
#

adata.obs['annotation'] = adata.obs['cell.type']

annota_default = Spectra.default_gene_sets.load()
annotations = {each:{} for each in Counter(adata.obs['annotation'])}
annotations['global'] = annota_default['global']
#

L = kst.estimate_L(adata, attribute = "annotation", highly_variable = True)
#run the model with default values
model = Spectra.est_spectra(
            adata=adata,
            gene_set_dictionary=annotations,
            use_highly_variable=True,
            cell_type_key="annotation",
            L = L,
            use_weights=True,
            lam=0.1, #varies depending on data and gene sets, try between 0.5 and 0.001
            delta=0.001,
            kappa=None,
            rho=0.001,
            use_cell_types=True,
            n_top_vals=50,
            label_factors=True,
            overlap_threshold=0.2,
            clean_gs = True,
            min_gs_num = 3,
            num_epochs=10000)
 
adata.write('MM.TME.235samples.TME.no_Eryth.Spectra.h5ad')
with open('MM.TME.235samples.TME.no_Eryth.Spectra.pickle','wb') as f:
        pickle.dump(model,f)

adata = sc.read_h5ad('MM.TME.235samples.TME.no_Eryth.Spectra.h5ad')
factors = adata.uns['SPECTRA_factors'] # factors x genes matrix that tells you how important each gene is to the resulting factors
markers = adata.uns['SPECTRA_markers'] # factors x n_top_vals list of n_top_vals top markers per factor
cell_scores = adata.obsm['SPECTRA_cell_scores'] # cells x factors matrix of cell scores
vocab = adata.var['spectra_vocab'] # boolean matrix of size # of genes that indicates the set of genes used to fit spectra
overlap = adata.uns['SPECTRA_overlap']
          
cell_scores = pd.DataFrame(cell_scores)
cell_scores.index = adata.obs_names
factors = pd.DataFrame(factors)
factors.columns = vocab.loc[vocab].index
factors.loc[0].nlargest(20)
#
cell_scores.columns = ['gene_module_'+str(each) for each in cell_scores.columns]
factors.index = ['gene_module_'+str(each) for each in factors.index]
markers = pd.DataFrame(markers)
markers.index = ['gene_module_'+str(each) for each in markers.index]

overlap = pd.DataFrame(overlap)
#
cell_scores.to_csv('MM.TME.235samples.TME.no_Eryth.Spectra.cell_scores.csv')
factors.to_csv('MM.TME.235samples.TME.no_Eryth.Spectra.gene_scores.csv')
markers.to_csv('MM.TME.235samples.TME.no_Eryth.Spectra.modified_marker_list.csv', header=False)
overlap.to_csv('MM.TME.235samples.TME.no_Eryth.Spectra.overlap.csv')

