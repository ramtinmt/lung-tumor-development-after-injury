"""scVI integration on CPU, with the sex genes removed before training.

Run from a terminal:   python -u integration_unbiased_sex.py
(or paste the whole file into one notebook cell).
"""
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'      # OpenMP conflict in this Anaconda install; must come before the imports

import numpy as np
import torch
import scanpy as sc
import scvi
import anndata as ad

FOLDER    = r'C:\Users\G20748795\Desktop\lab\Tumor'
INPUT     = 'preprocessed_sex.h5ad'
OUTPUT    = 'integrated20_sex.h5ad'
MODEL_DIR = 'scvi_n20_sex'
SEX_GENES = ['Xist', 'Tsix', 'Ddx3y', 'Eif2s3y', 'Kdm5d', 'Uty']

os.chdir(FOLDER)

# keep Windows from sleeping while this script runs (released automatically when it exits)
import ctypes
ctypes.windll.kernel32.SetThreadExecutionState(0x80000000 | 0x00000001)   # ES_CONTINUOUS | ES_SYSTEM_REQUIRED
torch.set_num_threads(os.cpu_count())             # all 24 cores
ad.settings.allow_write_nullable_strings = True   # otherwise the final write fails on this PC

adata = sc.read_h5ad(INPUT)
adata = adata[:, ~adata.var_names.isin(SEX_GENES)].copy()

scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key='sample')
model = scvi.model.SCVI(adata, n_latent=20)
model.train()

# save the small, valuable things first, so a problem with the big file cannot lose them
model.save(MODEL_DIR, overwrite=True)
latent = model.get_latent_representation()
np.save(os.path.join(MODEL_DIR, 'X_scVI.npy'), latent)

adata.obsm['X_scVI'] = latent
adata.write(OUTPUT)
