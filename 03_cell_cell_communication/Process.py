import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import json
from PIL import Image
import numpy as np

# === Load Visium data ===
adata = sc.read_visium(
    path="/gpfs/fs7/aafc/phenocart/PhenomicsProjects/RNASeq/SpatialTrans/Data/sample1a1",
    count_file="filtered_feature_bc_matrix.h5"
)

# === Preprocessing and Leiden clustering ===
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata)

# === Map clusters to biological cell types ===
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

# === Load scalefactors and image ===
with open("/gpfs/fs7/aafc/phenocart/PhenomicsProjects/RNASeq/SpatialTrans/Data/sample1a1/spatial/scalefactors_json.json") as f:
    scalefactors = json.load(f)

lowres_scale = scalefactors["tissue_lowres_scalef"]
img_path = "/gpfs/fs7/aafc/phenocart/PhenomicsProjects/RNASeq/SpatialTrans/Data/sample1a1/spatial/tissue_lowres_image.png"
img = Image.open(img_path)

# === Scale spatial coordinates ===
coords = pd.DataFrame(adata.obsm["spatial"], columns=["x", "y"], index=adata.obs.index)
coords["x"] = coords["x"] * lowres_scale
coords["y"] = coords["y"] * lowres_scale
cell_types = adata.obs["cell_type"]

# === Plot image with spot overlay ===
plt.figure(figsize=(10, 10))
plt.imshow(img, aspect='equal')
for ct in cell_types.unique():
    idx = cell_types == ct
    plt.scatter(
        coords.loc[idx, "x"],
        coords.loc[idx, "y"],
        s=100,            # Adjust spot size
        label=ct,
        alpha=0.25
    )

plt.gca().invert_yaxis()
plt.axis("off")
plt.legend(title="Cell Type", loc="upper right", bbox_to_anchor=(1.2, 1.0))
plt.title("Cell Type Overlay on H&E", fontsize=14)
plt.tight_layout()
plt.savefig("sample1a1_celltype_overlay_on_image.png", dpi=300, bbox_inches="tight")
plt.show()
