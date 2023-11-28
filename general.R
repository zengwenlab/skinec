seu <- readRDS("~/skin_project/seu_skin_SCT_harmony.rds")
library(Seurat)
library(ggplot2)
library(clusterProfiler)

table(seu$orig.ident)

DefaultAssay(seu)

seu <- RunUMAP(seu, dims = 1:50, verbose = T, reduction = "harmony")
seu <- FindNeighbors(seu, dims = 1:50, verbose = T, reduction = "harmony", assay = "SCT")
seu <- FindClusters(seu, verbose = T, resolution = seq(from = 0, by = 0.05, length = 10))

DimPlot(seu, raster = F, group.by = "SCT_snn_res.0.1", label = T, 
        repel = T, label.size = 8)

VlnPlot(seu, features = "PECAM1")
FeaturePlot(seu, features = "PECAM1")

library(clustree)

pdf(file = "clustree2.pdf", height = 10, width = 8)
clustree(seu) #+ coord_flip()
dev.off()


seu$seurat_clusters <- seu$SCT_snn_res.0.1
Idents(seu) <- seu$seurat_clusters



VlnPlot(seu, features = "LYVE1", pt.size = 0)

seu$seurat_clusters <- factor(seu$seurat_clusters, 
                              levels = c("0", "1", "2", "3", "4", "5", 
                                         "6", "7", "8", "9", "10", "11", 
                                         "12", "13", "14", "15", "16"))


Idents(seu)

degs_seu <- FindAllMarkers(seu, logfc.threshold = 0.25, min.pct = 0.5, only.pos = T)

library(dplyr)
top10 <- degs_seu %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
write.csv(top10, file="top10_cell_markers.csv")
write.csv(degs_seu, file = "degs_seu.csv")


marker_celltype <- c("COL1A1", "DCN", "KRTDAP", "KRT1", "KRT14", "KRT5", 
                     "CXCR4", "PTPRC", "ACTA2", "TAGLN", 
                     "SELE", "PECAM1", "LYZ", "HLA-DRA", 
                     "S100A8", "S100A9", "DCT", "TYRP1", "TPSAB1", "TPSB2", 
                     "CCL19", "CXCL12", "CCL21", "LYVE1", "PMP22", "CRYAB", 
                     "DES", "ACTG2", "COCH", "PTGDS", 
                     "HBB", "HBA2", "DCD", "SCGB1B2P")

marker_celltype2 <- c("KRTDAP", "KRT1", "KRT14", "KRT5", "COL1A1", 
                      "DCN",  "CCL19", "CXCL12", "COCH", 
                      "PTGDS", "ACTA2", "TAGLN", "DES", "ACTG2", 
                      "SELE", "PECAM1", "CCL21", "LYVE1", 
                      "CXCR4", "PTPRC", "S100A8", "S100A9", "LYZ", 
                      "HLA-DRA", "TPSAB1", "TPSB2", "DCT", 
                      "TYRP1", "PMP22", "CRYAB", "HBB", "HBA2", "DCD", "SCGB1B2P")

seu@meta.data$celltype <- seu@meta.data$seurat_clusters
Idents(seu) <- "celltype"
seu <- RenameIdents(seu, '0' = "Fib1", '1' = "Epi1", '2' = "Epi2", 
                    '3' = "Imm1", '4' = "Smc1", '5' = "VEc", 
                    '6' = "Mac", '7' = "Imm2", '8' = "Melan", 
                    '9' = "Mast", '10' = "Fib2", '11' = "LEc", 
                    '12' = "Swan", '13' = "Smc2", '14' = "Fib3", 
                    '15' = "Red", '16' = "Gland")
seu$celltype <- Idents(seu)
table(seu$celltype)

seu$celltype <- factor(seu$celltype, levels = c("Epi1", "Epi2", "Fib1", 
                                                "Fib2", "Fib3", "Smc1", 
                                                "Smc2", "VEc", "LEc", "Imm1", 
                                                "Imm2", "Mac","Mast", "Melan", 
                                                "Swan", "Red", "Gland"))


pdf(file = "celltype-marker.pdf", width = 13, height = 5)
DotPlot(seu, features = marker_celltype2, assay = "RNA", 
        group.by = "celltype", cols = c("gray90", "red")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "italic"))
dev.off()

saveRDS(seu, file = "seu_skin_all_711.rds")

library(randomcoloR)
palette <- randomColor(count = 60)  
palette <- distinctColorPalette(17) 
cols_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Paired")

pdf(file = "umap-celltype.pdf", height = 10, width = 12)
DimPlot(seu, reduction = "umap", group.by = "celltype", label = T, label.size = 8, repel = T, cols = palette)
dev.off()


pdf(file = "umap-Dataset.pdf", height = 10, width = 12)
DimPlot(seu, reduction = "umap", group.by = "Dataset", label = T, label.size = 7, repel = T, cols = palette)
dev.off()


pdf(file = "vln-cd31.pdf", width = 10, height = 7)
VlnPlot(seu, features = c("PECAM1", "LYVE1"), group.by = "celltype", pt.size = 0, ncol = 1)
dev.off()


rm(Degs_sig, Degs_sig_no_dup, top5);gc()


seu$ec.count <- ifelse(grepl("VEc", seu$celltype), "Ec", "nEc")

table(seu$ec.count)


proportion_barplot <- function(pbmc, cluster1, cluster2, col.cluster1){
  proportion <- as.data.frame(prop.table(x = table(pbmc@meta.data[[cluster1]], 
                                                   pbmc@meta.data[[cluster2]]), 
                                         margin = 2))
  colnames(proportion) <- c(cluster1, cluster2, "proportion")
  ggplot(proportion, aes_string(x = cluster2, group = cluster1)) + 
    geom_bar(aes_string(y="proportion", fill=cluster1), stat = "identity", width = 0.8)+ 
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) + 
    scale_y_continuous(limits=c(0, 1)) + 
    scale_x_discrete(limits = levels(pbmc@meta.data[[cluster2]])) + 
    scale_fill_manual(values = col.cluster1) +  
    cowplot::theme_cowplot() +
    ggtitle(label = "Proportion of cluster") # + coord_flip()
}

cols_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
cols_Set1 <- RColorBrewer::brewer.pal(n = 9, name = "Set1")

pdf(file = "ec_proportion.pdf", width = 20, height = 5)
proportion_barplot(seu, cluster1 = "ec.count", cluster2 = "orig.ident", 
                   col.cluster1 = cols_Set1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.3, vjust = 0.5))
dev.off()


pdf(file = "ec_proportion_gender.pdf", width = 5, height = 5)
proportion_barplot(seu, cluster1 = "ec.count", cluster2 = "Gender", 
                   col.cluster1 = cols_Set1)
dev.off()

seu$Age2 <- factor(seu$Age2, levels = c("Y", "M", "O"))
pdf(file = "ec_proportion_Age2.pdf", width = 5, height = 5)
proportion_barplot(seu, cluster1 = "ec.count", cluster2 = "Age2", 
                   col.cluster1 = cols_Set1)
dev.off()


pdf(file = "ec_proportion_location.pdf", width = 5, height = 5)
proportion_barplot(seu, cluster1 = "ec.count", cluster2 = "Location",
                   col.cluster1 = cols_Set1)
dev.off()


pdf(file = "ec_proportion_race.pdf", width = 5, height = 5)
proportion_barplot(seu, cluster1 = "ec.count", cluster2 = "Race", 
                   col.cluster1 = cols_Set1)
dev.off()

count_table <- table(seu$seurat_clusters, seu$orig.ident, seu$Dataset)
count_table <- as.data.frame(count_table)
colnames(count_table) <- c("seurat_clusters", "orig.ident", "Dataset", "Freq")

library(ggstatsplot)

pdf(file = "figS1-cluster-bar-harmony-dataset.pdf", width = 27, height = 8)
ggbarstats(count_table, x = Dataset, y = seurat_clusters, counts = Freq, package = "pals", palette = "polychrome", results.subtitle = F)
dev.off()


library(randomcoloR)
palette <- distinctColorPalette(65)
pdf(file = "figS1-cluster-bar-harmony-origident.pdf", width = 27, height = 8)
ggbarstats(count_table, x = orig.ident, y = seurat_clusters, counts = Freq, package = "khroma", palette = "stratigraphy", results.subtitle = F)
dev.off()



lisi_res_harmony <- lisi::compute_lisi(Embeddings(seu,reduction = "umap"), 
                                       seu@meta.data, 
                                       c('seurat_clusters', 
                                         'orig.ident', 'Dataset'))

p61 <- ggplot(lisi_res_harmony, aes(x=orig.ident))+
  geom_histogram(aes(y=after_stat(density)),
                 binwidth = .5,
                 colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  geom_vline(aes(xintercept=mean(orig.ident, na.rm=T)),
             color="red", linetype="dashed", linewidth=1)

p62 <- ggplot(lisi_res_harmony, aes(x=Dataset))+
  geom_histogram(aes(y=after_stat(density)),
                 binwidth = .5,
                 colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  geom_vline(aes(xintercept=mean(Dataset, na.rm=T)),
             color="red", linetype="dashed", linewidth=1)

p64 <- ggplot(lisi_res_harmony, aes(x=seurat_clusters))+
  geom_histogram(aes(y=after_stat(density)),
                 binwidth = .5,
                 colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  geom_vline(aes(xintercept=mean(seurat_clusters, na.rm=T)),
             color="red", linetype="dashed", linewidth=1) +
  scale_x_continuous(breaks = scales::breaks_width(1))


p6 <- ((p61|p62|p64)) + theme_classic()
ggsave(p6, filename = "figS1-score-harmony.pdf", width = 15, height = 4)


library(randomcoloR)
palette <- distinctColorPalette(65)
pdf(file = "umap-orig.idents.pdf", width = 13, height = 8)
DimPlot(seu, reduction = "umap", label = T, repel = T, group.by = "orig.ident", 
        cols = palette)
dev.off()


seu <- RunUMAP(seu, dims = 1:50, verbose = T, reduction = "pca")
seu <- FindNeighbors(seu, dims = 1:50, verbose = T, reduction = "pca")
seu <- FindClusters(seu, verbose = T, resolution = 0.1)

count_table <- table(seu$seurat_clusters, seu$orig.ident, seu$Dataset)
count_table <- as.data.frame(count_table)
colnames(count_table) <- c("seurat_clusters", "orig.ident", "Dataset", "Freq")

library(ggstatsplot)

pdf(file = "figS1-cluster-bar-unharmony-dataset.pdf", width = 27, height = 8)
ggbarstats(count_table, x = Dataset, y = seurat_clusters, counts = Freq, 
           package = "pals", palette = "polychrome", results.subtitle = F)
dev.off()

library(randomcoloR)
palette <- distinctColorPalette(65)
pdf(file = "figS1-cluster-bar-unharmony-origident.pdf", width = 27, height = 8)
ggbarstats(count_table, x = orig.ident, y = seurat_clusters, counts = Freq, package = "khroma", palette = "stratigraphy", results.subtitle = F)
dev.off()



lisi_res_unharmony <- lisi::compute_lisi(Embeddings(seu,reduction = "umap"), 
                                         seu@meta.data, 
                                         c('seurat_clusters', 
                                           'orig.ident', 'Dataset'))

p61 <- ggplot(lisi_res_unharmony, aes(x=orig.ident))+
  geom_histogram(aes(y=after_stat(density)),
                 binwidth = .5,
                 colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  geom_vline(aes(xintercept=mean(orig.ident, na.rm=T)),
             color="red", linetype="dashed", linewidth=1)

p62 <- ggplot(lisi_res_unharmony, aes(x=Dataset))+
  geom_histogram(aes(y=after_stat(density)),
                 binwidth = .5,
                 colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  geom_vline(aes(xintercept=mean(Dataset, na.rm=T)),
             color="red", linetype="dashed", linewidth=1)

p64 <- ggplot(lisi_res_unharmony, aes(x=seurat_clusters))+
  geom_histogram(aes(y=after_stat(density)),
                 binwidth = .5,
                 colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  geom_vline(aes(xintercept=mean(seurat_clusters, na.rm=T)),
             color="red", linetype="dashed", linewidth=1) +
  scale_x_continuous(breaks = scales::breaks_width(1))


p6 <- ((p61|p62|p64)) + theme_classic()
ggsave(p6, filename = "figS1-score-unharmony.pdf", width = 15, height = 4)


pdf(file = "umap-gender.pdf", height = 10, width = 13)
DimPlot(seu, group.by = "Gender", reduction = "umap", label = T, repel = T)
dev.off()


pdf(file = "umap-age.pdf", height = 10, width = 12)
DimPlot(seu, group.by = "Age2", reduction = "umap", label = T, repel = T)
dev.off()


pdf(file = "umap-location.pdf", height = 10, width = 12)
DimPlot(seu, group.by = "Location", reduction = "umap", label = T, repel = T)
dev.off()


pdf(file = "umap-race.pdf", height = 10, width = 12)
DimPlot(seu, group.by = "Race", reduction = "umap", label = T, repel = T)
dev.off()


count_table <- table(seu$seurat_clusters, seu$Gender, seu$Age2, seu$Location, seu$Race)
count_table <- as.data.frame(count_table)
colnames(count_table) <- c("seurat_clusters", "Gender", "Age2", "Location", "Race", "Freq")

library(ggstatsplot)

pdf(file = "figS1-cluster-bar-harmony-Gender.pdf", width = 15, height = 8)
ggbarstats(count_table, x = Gender, y = seurat_clusters, counts = Freq, package = "pals", palette = "polychrome", results.subtitle = F)
dev.off()


pdf(file = "figS1-cluster-bar-harmony-Age.pdf", width = 15, height = 8)
ggbarstats(count_table, x = Age2, y = seurat_clusters, counts = Freq, package = "pals", palette = "polychrome", results.subtitle = F)
dev.off()


pdf(file = "figS1-cluster-bar-harmony-Location.pdf", width = 15, height = 8)
ggbarstats(count_table, x = Location, y = seurat_clusters, counts = Freq, package = "pals", palette = "polychrome", results.subtitle = F)
dev.off()


pdf(file = "figS1-cluster-bar-harmony-Race.pdf", width = 15, height = 8)
ggbarstats(count_table, x = Race, y = seurat_clusters, counts = Freq, package = "pals", palette = "polychrome", results.subtitle = F)
dev.off()




