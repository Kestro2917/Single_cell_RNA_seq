# installation of "Seurat" package
install.packages("Seurat")

# installation of "SeuratDisk" package
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

# Load "Seurat" and "SeuratDisk" libraries

library("Seurat")
library("SeuratDisk")

# 10x CellRanger .HDF5 format

hdf5_obj <- Read10X_h5(filename = "C:\\Users\\biksk\\OneDrive\\Desktop\\Desktop\\Single_cell_RNA_seq\\Data\\20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
hdf5_obj[1:10, 1:10]

seurat_hdf5 <- CreateSeuratObject(counts = hdf5_obj)
seurat_hdf5
str(seurat_hdf5)

# .mtx file

mtx_obj <- ReadMtx(mtx = "C:\\Users\\biksk\\OneDrive\\Desktop\\Desktop\\Single_cell_RNA_seq\\Data\\20k_PBMC_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix\\raw_feature_bc_matrix\\matrix.mtx", features = "C:\\Users\\biksk\\OneDrive\\Desktop\\Desktop\\Single_cell_RNA_seq\\Data\\20k_PBMC_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix\\raw_feature_bc_matrix\\features.tsv", cells = "C:\\Users\\biksk\\OneDrive\\Desktop\\Desktop\\Single_cell_RNA_seq\\Data\\20k_PBMC_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix\\raw_feature_bc_matrix\\barcodes.tsv")


# .h5ad formate
# Step 1: convert AnnData object to an h5Seurat file

Convert("C:\\Users\\biksk\\OneDrive\\Desktop\\Desktop\\Single_cell_RNA_seq\\Data\\adata_SS2_for_download.h5ad", dest = "h5seurat", overwrite = TRUE)

# Step 2: Load h5Seurat file into a Seurat object
seurate_anndata <- LoadH5Seurat("C:\\Users\\biksk\\OneDrive\\Desktop\\Desktop\\Single_cell_RNA_seq\\Data\\adata_SS2_for_download.h5seurat")
str(seurate_anndata)