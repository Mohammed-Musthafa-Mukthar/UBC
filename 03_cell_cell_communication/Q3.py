import scanpy as sc
import pandas as pd
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
from ktplotspy import plot_cpdb_heatmap

#1: Load Visium Data 
adata = sc.read_visium(
    path="/gpfs/fs7/aafc/phenocart/PhenomicsProjects/RNASeq/SpatialTrans/Data/sample1a1",
    count_file="filtered_feature_bc_matrix.h5"
)

# 2: Preprocessing and Clustering
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata)

# 3: Annotate Cell Types
cluster_to_type = {
    '0': 'SMC',
    '1': 'Endothelial',
    '2': 'Fibroblast',
    '3': 'Adipocyte',
    '4': 'SMC',
    '5': 'Macrophage',
    '6': 'Fibroblast'
}
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_to_type)

# 4: Export meta.txt
meta = pd.DataFrame({
    'Cell': adata.obs_names,
    'cell_type': adata.obs['cell_type']
})
meta_file_path = "meta.txt"
meta.to_csv(meta_file_path, sep='\t', index=False)

# 5: Export counts.txt
expr_array = adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X
expr = pd.DataFrame(
    expr_array.T,  # Transpose to match CellPhoneDB format
    index=adata.var_names,
    columns=adata.obs_names
)
expr.insert(0, "Gene", expr.index)  # Add Gene name column
counts_file_path = "counts.txt"
expr.to_csv(counts_file_path, sep='\t', index=False)

print("meta.txt and counts.txt generated successfully.")

# 6: Run CellPhoneDB
cpdb_file_path = "cellphonedb.zip"  # Provide the correct path if needed
output_path = "out"
os.makedirs(output_path, exist_ok=True)

print("Running CellPhoneDB...")
cpdb_statistical_analysis_method.call(
    cpdb_file_path=cpdb_file_path,
    meta_file_path=meta_file_path,
    counts_file_path=counts_file_path,
    counts_data='hgnc_symbol',
    output_path=output_path,
    threshold=0.1,
    iterations=1000,
    threads=4,
    debug_seed=-1,
    result_precision=3
)
print("CellPhoneDB analysis complete.")

# 7: Visualize Significant Interactions
print("Creating heatmap...")

pvalues_file = sorted(glob.glob(os.path.join(output_path, "statistical_analysis_pvalues_*.txt")))[-1]
pvals = pd.read_csv(pvalues_file, sep="\t")

fig = plot_cpdb_heatmap(
    pvals=pvals,
    degs_analysis=False,
    figsize=(10, 10),
    title="Sum of significant interactions"
)
heatmap_path = os.path.join(output_path, "cellphonedb_heatmap.png")
fig.savefig(heatmap_path, dpi=300)
print("Heatmap saved to:", heatmap_path)
