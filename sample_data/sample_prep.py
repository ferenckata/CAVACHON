''' Custom script to do stratified filting of the cells in the example data
so that it can be used for quick demonstration'''

import scanpy
import math

# read in data
atac_mod = scanpy.read_h5ad("Chen-2019-ATAC.h5ad")
rna_mod = scanpy.read_h5ad("Chen-2019-RNA.h5ad")


def filter_cells(
        atac_adata: scanpy.AnnData,
        rna_adata: scanpy.AnnData,
        target_percentage: float = 0.1) -> scanpy.AnnData:
    ''' Downsample each cell type to x% of its original size '''
    cluster_key = "cell_type"
    grouped = atac_adata.obs.groupby(cluster_key)
    downsampled_indices = []

    for _, group in grouped:
        down_size = math.ceil(group.size * target_percentage)
        downsampled_indices.extend(group.sample(down_size).index)

    atac_downsampled = atac_adata[downsampled_indices]
    rna_downsampled = rna_adata[downsampled_indices]
    return atac_downsampled, rna_downsampled


prop_cells_per_type = 0.1
downsampled_atac, downsampled_rna = filter_cells(
    atac_adata=atac_mod,
    rna_adata=rna_mod,
    target_percentage=prop_cells_per_type)

downsampled_atac.write_h5ad(filename="Chen-2019-ATAC_10pc.h5ad")
downsampled_rna.write_h5ad(filename="Chen-2019-RNA_10pc.h5ad")
