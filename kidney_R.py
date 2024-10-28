import os

import anndata
import scanpy as sc
from matplotlib import rcParams
from networkx.algorithms.bipartite import biadjacency_matrix

import scglue


import numpy as np
import scipy.io
# %%
scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)

PATH = "mJSNMF"
os.makedirs(PATH, exist_ok=True)

# %% [markdown]
# # Read data

# %%
mat = scipy.io.loadmat('./mJSNMF/kidney_5k_10k.mat')
X1 = mat['X1']
X2 = mat['X2']
label = mat['label'][:,0]
import pandas as pd
labels=[pd.DataFrame(label)[0][i][0] for i in range(len(label))]

# create AnnData object
rna = anndata.AnnData(X1.T,dtype = np.float64)
atac = anndata.AnnData(X2.T,dtype = np.float64)
rna.obs['celltype'] = label

genes = np.array([ val[0]  for val in mat['genes'][:,0]])
peaks = np.array([ val[0]  for val in mat['peaks'][:,0]])
rna.var_names=genes
scglue.data.get_gene_annotation(
    rna, gtf="gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz",
    gtf_by="gene_name"
)
atac.var_names=peaks
split = atac.var_names.str.split(r"[:-]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[-2]).astype(int)
atac.var["chromEnd"] = split.map(lambda x: x[-1]).astype(int)
# %%
# unique_labels = list(set(rna.var["chrom"]))
used_chroms = {f"chr{x}" for x in range(1, 20)}.union({"chrX"}).union({"chrY"})
#
# %%
rna = rna[:,[item in used_chroms for item in rna.var["chrom"]]]
# sc.pp.filter_genes(rna, min_counts=1)
# rna.obs_names += "-RNA"

# %%
atac = atac[:,[item in used_chroms for item in atac.var["chrom"]]]
# sc.pp.filter_genes(atac, min_counts=1)
# atac.obs_names += "-ATAC"


genes = scglue.genomics.Bed(rna.var.assign(name=rna.var_names))
peaks = scglue.genomics.Bed(atac.var.assign(name=atac.var_names))
overlap_graph = scglue.genomics.window_graph(
    genes.expand(2000, 0), peaks, 0,
    attr_fn=lambda l, r, d: {
        "weight": 1.0,
        "type": "overlap"
    }
)

A=biadjacency_matrix(overlap_graph,genes.index,peaks.index)
A.nnz

from scipy import sparse
sparse.save_npz('./mJSNMF/kidney_R.npz', A)  # 保存
rna.write("./mJSNMF/kidney_5k_10k_rna.h5ad")
atac.write("./mJSNMF/kidney_5k_10k_atac.h5ad")
