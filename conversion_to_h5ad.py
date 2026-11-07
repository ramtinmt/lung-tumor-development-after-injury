import scanpy as sc
from pathlib import Path

path = Path('/Users/G20748795/Desktop/lab/Tumor/')

for file in path.glob('*.h5'):
    adata = sc.read_10x_h5(file)
    adata.write_h5ad(file.with_suffix('.h5ad'))
print("All h5 files converted to h5ad")


for file in path.glob('*.loom'):
    adata = sc.read_loom(file)
    adata.write_h5ad(file.with_suffix('.h5ad'))
print("All loom files converted to h5ad")
