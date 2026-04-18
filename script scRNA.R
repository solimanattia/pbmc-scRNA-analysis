# ================================================================
# Title: scRNA-seq Analysis of PBMC 3k Dataset
# Author: soliman attia 
# Date: 2026-04-17
# Description: Full Seurat pipeline from raw data to cell annotation and DEGs
# ================================================================

# 1. تحميل الباكدجات -----------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# 2. تحميل الداتا -------------------------------------------------
# فك الضغط من pbmc3k_filtered_gene_bc_matrices.tar.gz الأول
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# 3. Quality Control ----------------------------------------------
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 4. Normalization + اختيار الجينات المتغيرة ----------------------
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# 5. Scaling + PCA ------------------------------------------------
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# 6. Clustering + UMAP --------------------------------------------
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

# 7. Cell Type Annotation -----------------------------------------
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", 
                     "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# 8. حفظ النتايج --------------------------------------------------
# UMAP ملون
png("umap_labeled.png", width=2000, height=1800, res=300)
DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 6) + NoLegend()


# حفظ الـ Seurat object
saveRDS(pbmc, file = "pbmc_annotated.rds")

# 9. Differential Expression --------------------------------------
deg_mono <- FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono")
write.csv(deg_mono, file = "DEG_CD14_vs_FCGR3A_mono.csv")

# 10. معلومات الجلسة عشان اللي يكرر التحليل -----------------------
sessionInfo()
