seu <- readRDS("~/skin_project/Figure1_all/seu_skin_all_711.rds")

seu_ec <- subset(seu, celltype == "VEc")

table(seu_ec$Gender)

seu_ec

DefaultAssay(seu_ec) <- "RNA"

seu_ec$ec.count <- NULL
seu_ec$celltype <- NULL
seu_ec$seurat_clusters <- NULL
seu_ec$SCT_snn_res.0.45 <- NULL
seu_ec$SCT_snn_res.0.4 <- NULL
seu_ec$SCT_snn_res.0.35 <- NULL
seu_ec$SCT_snn_res.0.3 <- NULL
seu_ec$SCT_snn_res.0.25 <- NULL
seu_ec$SCT_snn_res.0.2 <- NULL
seu_ec$SCT_snn_res.0.15 <- NULL
seu_ec$SCT_snn_res.0.1 <- NULL
seu_ec$SCT_snn_res.0.05 <- NULL
seu_ec$SCT_snn_res.0 <- NULL
seu_ec$percent.mt <- NULL

seu_ec$S.Score <- NULL
seu_ec$G2M.Score <- NULL
seu_ec$Phase <- NULL
seu_ec@reductions$pca <- NULL
seu_ec@reductions$harmony <- NULL
seu_ec@reductions$umap <- NULL


###
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
seu_ec <- CellCycleScoring(seu_ec, s.features = s.genes, g2m.features = g2m.genes)
seu_ec <- PercentageFeatureSet(seu_ec, pattern = "^MT-", col.name = "percent.mt")
seu_ec <- SCTransform(seu_ec, 
                      vars.to.regress = c("nCount_RNA", "percent.mt", 
                                          "S.Score", "G2M.Score"), verbose = T)

DefaultAssay(seu_ec)
seu <- RunPCA(seu_ec, verbose = T)

library(harmony)
seu_ec <- RunHarmony(seu_ec, group.by.vars = "orig.ident", max.iter.harmony = 20)

saveRDS(seu_ec, file = "seu_ec_SCT_harmony.rds")



###
seu_ec <- RunUMAP(seu_ec, dims = 1:50, verbose = T, reduction = "harmony")
seu_ec <- FindNeighbors(seu_ec, dims = 1:50, verbose = T, reduction = "harmony")
seu_ec <- FindClusters(seu_ec, verbose = T, 
                       resolution = seq(from = 0, by = 0.05, length = 10))


####

library(clustree)

pdf(file = "clustree2.pdf", height = 10, width = 8)
clustree(seu_ec) #+ coord_flip()
dev.off()


###

DimPlot(seu_ec, raster = F, reduction = "umap", label = T, group.by = "SCT_snn_res.0.15")


###

table(seu_ec$SCT_snn_res.0.15)

seu_ec$SCT_snn_res.0.15 <- factor(seu_ec$SCT_snn_res.0.15, 
                                  levels = c("0", "1", "2", "3", "4", "5", 
                                             "6", "7", "8", "9", "10"))
seu_ec$seurat_clusters <- seu_ec$SCT_snn_res.0.15
Idents(seu_ec) <- seu_ec$seurat_clusters

pdf(file = "umap.pdf", width = 7, height = 5)
DimPlot(seu_ec, raster = F, reduction = "umap", label = T, group.by = "seurat_clusters")
dev.off()


###
saveRDS(seu_ec, file = "seu_ec_715.rds")

####
Idents(seu_ec) <- "seurat_clusters"
degs_ec <- FindAllMarkers(seu_ec, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)
write.csv(degs_ec, file = "degs_ec.csv")

###

library(dplyr)
Degs_sig <- subset(degs_ec, p_val_adj < 0.01)
Degs_sig %>% group_by(gene) %>% top_n(1, avg_log2FC) -> Degs_sig_no_dup
Degs_sig_no_dup  %>% group_by(cluster) %>% top_n(3, -p_val_adj) -> top3
top3 %>% group_by(cluster) %>% top_n(3, avg_log2FC) -> top3

library(scRNAtoolVis)
cols_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
pdf(file = "heatmap.ec.pdf", width = 8, height = 9)
AverageHeatmap(object = seu_ec, group.by = "seurat_clusters", annoCol = T, 
               myanCol = cols_clusters[2:12],
               markerGene = Degs_sig_no_dup$gene, row_title = "", 
               showRowNames = F, markGenes = top3$gene)
dev.off()



###

degs_ec2 <- FindAllMarkers(seu_ec, min.pct = 0.3, logfc.threshold = 0.5, only.pos = T)
Degs_sig_subset <- subset(degs_ec2, p_val_adj < 0.01)
enrich_data <- Degs_sig_subset[, c("cluster","gene")]
library(org.Hs.eg.db)
library(clusterProfiler)
go_subset_id <- bitr(enrich_data[,2],
                     fromType="SYMBOL",
                     toType=c("ENTREZID","ENSEMBL","SYMBOL"),
                     OrgDb="org.Hs.eg.db")
library(dplyr)
go_subset <- inner_join(enrich_data,go_subset_id,by=c("gene"="SYMBOL"))
GO_CompareCluster_Reslut <- compareCluster(gene~cluster, 
                                           data=go_subset, 
                                           fun="enrichGO", 
                                           OrgDb = org.Hs.eg.db, 
                                           keyType = 'SYMBOL', 
                                           ont = "BP", 
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05, 
                                           qvalueCutoff  = 0.2)

GO_filtered <- clusterProfiler::simplify(GO_CompareCluster_Reslut, 
                                         cutoff = 0.5, 
                                         by = "p.adjust")
GO_filtered@compareClusterResult$cluster <- factor(GO_filtered@compareClusterResult$cluster, 
                                                   levels = c("0", "1", "2", "3", 
                                                              "4", "5", "6", 
                                                              "7", "8", "9", "10"))
pdf(file = "go-ec-bp2.pdf", height = 12, width = 10)
dotplot(GO_filtered, x = "cluster", color = "p.adjust", showCategory = 6, 
        label_format=150)
dev.off()


##

GO_CompareCluster_Reslut <- compareCluster(gene~cluster, 
                                           data=go_subset, 
                                           fun="enrichGO", 
                                           OrgDb = org.Hs.eg.db, 
                                           keyType = 'SYMBOL', 
                                           ont = "MF", 
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05, 
                                           qvalueCutoff  = 0.2)
GO_filtered <- clusterProfiler::simplify(GO_CompareCluster_Reslut, 
                                         cutoff = 0.5, 
                                         by = "p.adjust")

pdf(file = "go-ec-mf2.pdf", height = 12, width = 10)
dotplot(GO_filtered, x = "cluster", color = "p.adjust", showCategory = 6, 
        label_format=150)
dev.off()


##
GO_CompareCluster_Reslut <- compareCluster(gene~cluster, 
                                           data=go_subset, 
                                           fun="enrichGO", 
                                           OrgDb = org.Hs.eg.db, 
                                           keyType = 'SYMBOL', 
                                           ont = "CC", 
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05, 
                                           qvalueCutoff  = 0.2)
GO_filtered <- clusterProfiler::simplify(GO_CompareCluster_Reslut, 
                                         cutoff = 0.5, 
                                         by = "p.adjust")

pdf(file = "go-ec-cc2.pdf", height = 12, width = 10)
dotplot(GO_filtered, x = "cluster", color = "p.adjust", showCategory = 6, label_format=150)
dev.off()

##
gene_arterial <- c("Dll4", "Igfbp3", "Unc5b", "Gja4", 
                   "Hey1", "Mecom", "Efnb2", "Epas1", 
                   "Vegfc", "Cxcr4", "Bmx",
                   "Nrp1", "Gja5")
gene_venous <- c("Nr2f2", "Nrp2","Aplnr")
gene_arterial <- toupper(gene_arterial)
gene_venous <- toupper(gene_venous)

seu_ec <- AddModuleScore(seu_ec, features = list(gene_arterial, gene_venous), 
                         name = c("Arterial", "Venous"))
colnames(seu_ec@meta.data)[28:29] <- c("Arterial_score", "Venous_score")

ggdata <- FetchData(seu_ec, vars = c("seurat_clusters", 
                                     "Arterial_score", "Venous_score"))
ggdata <- reshape2::melt(ggdata, id = "seurat_clusters")
colnames(ggdata) <- c("seurat_clusters","feature","expression")
cols_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
library(cowplot)

pdf(file = "A_V_score.pdf", width = 10, height = 4)

ggplot(ggdata, aes(seurat_clusters, expression, fill = seurat_clusters)) + 
  geom_violin(scale = "width")+ 
  theme_cowplot()+ 
  ggtitle(label = "Ec characteristics and signalings")+
  geom_boxplot(outlier.size = 0.5, notch = F, fill = "white", width = 0.1) + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank())+
  facet_wrap(~feature, scales = "free_y", ncol = 4)  + 
  scale_fill_manual(values = cols_clusters[2:12])

dev.off()


##
pdf(file = "maker.pdf", width = 8, height = 10)
DotPlot(seu_ec, features = c("SEMA3G", "PLVAP", "SELE", "ACKR1", "FBLN2"))
dev.off()









