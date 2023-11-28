seu_ec <- readRDS("~/skin_project/figure2ec/seu_ec_715.rds")

library(Seurat)
library(ggplot2)

table(seu_ec$Location)

seu_ec <- subset(seu_ec, Location %in% c("abdomen", "arms&legs", "chest", "head"))
seu_ec$orig.ident <- droplevels(seu_ec$orig.ident)
table(seu_ec$orig.ident)


pdf(file = "umap-loc.pdf", height = 5, width = 10)
DimPlot(seu_ec, group.by = "seurat_clusters", split.by = "Location", label = T, repel = T)
dev.off()


proportion_barplot <- function(pbmc, cluster1, cluster2, col.cluster1){
  proportion <- as.data.frame(prop.table(x = table(pbmc@meta.data[[cluster1]], pbmc@meta.data[[cluster2]]), margin = 2))
  colnames(proportion) <- c(cluster1, cluster2, "proportion")
  ggplot(proportion, aes_string(x = cluster2, group = cluster1)) + 
    geom_bar(aes_string(y="proportion", fill=cluster1), stat = "identity", width = 0.8)+ 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_y_continuous(limits=c(0, 1)) + 
    scale_x_discrete(limits = levels(pbmc@meta.data[[cluster2]])) + 
    scale_fill_manual(values = col.cluster1) +  
    cowplot::theme_cowplot() +
    ggtitle(label = "Proportion of cluster") # + coord_flip()
}

cols_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
pdf(file = "Location-propotion.pdf", width = 8, height = 9)
proportion_barplot(seu_ec, cluster1 = "seurat_clusters", 
                   cluster2 = "Location", col.cluster1 = cols_clusters)
dev.off()


Idents(seu_ec) <- "Location"

seu_ec_deg_location <- FindAllMarkers(seu_ec, only.pos = T)

Degs_sig <- subset(seu_ec_deg_location, p_val_adj < 0.01)
library(dplyr)
Degs_sig %>% group_by(gene) %>% top_n(5, avg_log2FC) -> Degs_sig_no_dup
Degs_sig_no_dup  %>% group_by(cluster) %>% top_n(5, -p_val_adj) -> top5
top5 %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5
cols_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
library(scRNAtoolVis)
pdf(file = "heatmap.ec-location.pdf", width = 7, height = 8)
AverageHeatmap(object = seu_ec, group.by = "Location", annoCol = T, myanCol = cols_clusters[2:5],
               markerGene = Degs_sig_no_dup$gene, row_title = "", 
               showRowNames = F, markGenes = top5$gene)
dev.off()


Degs_sig <- subset(seu_ec_deg_location, p_val_adj < 0.01)
write_csv(Degs_sig, file = "Degs_sig.csv")

enrich_data <- Degs_sig[, c("cluster","gene")]
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
                                           ont = "ALL", 
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05, 
                                           qvalueCutoff  = 0.2)
GO_filtered <- simplify(GO_CompareCluster_Reslut, 
                        cutoff = 0.5, 
                        by = "p.adjust")

pdf(file = "go-ec-location-all.pdf", height = 15, width = 12)
dotplot(GO_filtered, x = "cluster", color = "p.adjust", 
        showCategory = 6, split = "ONTOLOGY", label_format=150) + 
  facet_grid(ONTOLOGY~., scales = "free")
dev.off()



GO_CompareCluster_Reslut <- compareCluster(gene~cluster, 
                                           data=go_subset, 
                                           fun="enrichGO", 
                                           OrgDb = org.Hs.eg.db, 
                                           keyType = 'SYMBOL', 
                                           ont = "BP", 
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05, 
                                           qvalueCutoff  = 0.2)
GO_filtered <- simplify(GO_CompareCluster_Reslut, 
                        cutoff = 0.5, 
                        by = "p.adjust")

pdf(file = "go-ec-location-bp.pdf", height = 10, width = 11)
dotplot(GO_filtered, x = "cluster", color = "p.adjust", showCategory = 10, label_format=90) 
dev.off()

pdf(file = "ec-location-marker.pdf", width = 8, height = 10)
VlnPlot(seu_ec, features = c("SEPTIN7", "RBIS", "CXCL2", "CCL2", 
                             "FMNL2", "NAV2", "CD200", "S100A2"), ncol = 2, pt.size = 0)
dev.off()


library(Seurat)
library(ggplot2)
library(irGSEA)
seu_ec <- subset(seu_ec, Location %in% c("abdomen", "arms&legs", "chest", "head"))
seu_ec$orig.ident <- droplevels(seu_ec$orig.ident)
Idents(seu_ec) <- "Location"

seu_ec <- irGSEA.score(object = seu_ec, assay = "SCT", 
                       slot = "data", seeds = 123, ncores = 10,
                       custom = F, geneset = NULL, msigdb = T, 
                       species = "Homo sapiens", category = "H",  
                       subcategory = NULL, geneid = "symbol",
                       method = c("AUCell", "UCell", "singscore", 
                                  "ssgsea", "JASMINE", "viper"),
                       kcdf = 'Gaussian')
saveRDS(seu_ec, file = "seu_ec_score_location_H.rds")
result.dge.H <- irGSEA.integrate(object = seu_ec, 
                                 group.by = "Location",
                                 metadata = NULL, col.name = NULL,
                                 method = c("AUCell","UCell","singscore",
                                            "ssgsea", "JASMINE", "viper"))
saveRDS(result.dge.H, file = "result.dge.location.H.rds")

####
print(msigdbr::msigdbr_collections(), n=30) 
seu_ec <- irGSEA.score(object = seu_ec, assay = "SCT", 
                       slot = "data", seeds = 123, ncores = 10,
                       custom = F, geneset = NULL, msigdb = T, 
                       species = "Homo sapiens", category = "C2",  
                       subcategory = "CP:KEGG", geneid = "symbol",
                       method = c("AUCell", "UCell", "singscore", 
                                  "ssgsea", "JASMINE", "viper"),
                       kcdf = 'Gaussian')
saveRDS(seu_ec, file = "seu_ec_score_location.KEGG.rds")
result.dge.KEGG <- irGSEA.integrate(object = seu_ec, 
                                    group.by = "Location",
                                    metadata = NULL, col.name = NULL,
                                    method = c("AUCell","UCell","singscore",
                                               "ssgsea", "JASMINE", "viper"))
saveRDS(result.dge.KEGG, file = "result.dge.location.KEGG.rds")


pdf(file = "heatmap-location-kegg-RRA.pdf", width = 15, height = 15)
irGSEA.heatmap(object = result.dge.location.KEGG, method = "RRA", top = 50, show.geneset = NULL, heatmap.width = 20, heatmap.heigh = 18)
dev.off()


pdf(file = "location-H-upset-RRA.pdf", width = 7, height = 5)
irGSEA.upset(object = result.dge.location.H, method = "RRA")
dev.off()



pdf(file = "location-kegg-upset-RRA.pdf", width = 7, height = 5)
irGSEA.upset(object = result.dge.location.KEGG, method = "RRA")
dev.off()

pdf(file = "location-kegg-barplot-H.pdf", height = 7, width = 13)
irGSEA.barplot(object = result.dge.location.H, method = c("AUCell", "UCell", "singscore", "ssgsea", "JASMINE", "viper", "RRA"))
dev.off()


pdf(file = "location-kegg-barplot-KEGG.pdf", height = 7, width = 13)
irGSEA.barplot(object = result.dge.location.KEGG, method = c("AUCell", "UCell", "singscore", "ssgsea", "JASMINE", "viper", "RRA"))
dev.off()

pdf(file = "location-halfvlnplot-NOTCH-aucell-H.pdf", width = 10, height = 7)
irGSEA.halfvlnplot(object = seu_ec_score_location_H, method = "AUCell", show.geneset = "HALLMARK-NOTCH-SIGNALING")
dev.off()

pdf(file = "location-halfvlnplot-ANGIOGENESIS-aucell-H.pdf", width = 10, height = 7)
irGSEA.halfvlnplot(object = seu_ec_score_location_H, method = "AUCell", show.geneset = "HALLMARK-ANGIOGENESIS")
dev.off()


pdf(file = "location-densityheatmap-aucell-VEGF.pdf", width = 10, height = 8)
irGSEA.densityheatmap(object = seu_ec_score_location.KEGG, method = "AUCell", show.geneset = "KEGG-VEGF-SIGNALING-PATHWAY")
dev.off()




