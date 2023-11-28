seu_ec <- readRDS("~/skin_project/figure2ec/seu_ec_715.rds")
library(Seurat)
library(ggplot2)
table(seu_ec$Race)

seu_ec_race <- subset(seu_ec, Race %in% c("African", "Asian", "Caucasian"))
seu_ec_race$orig.ident <- droplevels(seu_ec_race$orig.ident)
seu_ec_race$Race <- droplevels(seu_ec_race$Race)


pdf(file = "umap-race.pdf", height = 5, width = 10)
DimPlot(seu_ec_race, group.by = "seurat_clusters", split.by = "Race", label = T, repel = T)
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
pdf(file = "Race-propotion.pdf", width = 8, height = 9)
proportion_barplot(seu_ec_race, cluster1 = "seurat_clusters", 
                   cluster2 = "Race", col.cluster1 = cols_clusters)
dev.off()

Idents(seu_ec_race) <- "Race"

seu_ec_deg_race <- FindAllMarkers(seu_ec_race, only.pos = T)

Degs_sig <- subset(seu_ec_deg_race, p_val_adj < 0.01)
library(dplyr)
Degs_sig %>% group_by(gene) %>% top_n(5, avg_log2FC) -> Degs_sig_no_dup
Degs_sig_no_dup  %>% group_by(cluster) %>% top_n(5, -p_val_adj) -> top5
top5 %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5
cols_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Paired")

seu_ec_race$Race <- factor(seu_ec_race$Race, levels = c("African", "Asian", "Caucasian"))
library(scRNAtoolVis)
pdf(file = "heatmap.ec-race.pdf", width = 7, height = 9)
AverageHeatmap(object = seu_ec_race, group.by = "Race", annoCol = T, myanCol = cols_clusters[2:4],
               markerGene = Degs_sig_no_dup$gene, row_title = "", 
               showRowNames = F, markGenes = top5$gene)
dev.off()

write_csv(Degs_sig, file = "Degs_sig_race.csv")

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

pdf(file = "go-ec-race-all.pdf", height = 13, width = 11)
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

pdf(file = "go-ec-race-bp.pdf", height = 9, width = 8)
dotplot(GO_filtered, x = "cluster", color = "p.adjust", showCategory = 10, label_format=90) 
dev.off()

pdf(file = "race-marker.pdf", height = 8, width = 10)
VlnPlot(seu_ec_race, features = c("RPL7", "RPL31", "HMOX1", "NAV2", "FABP4", "MKL2"), pt.size = 0)
dev.off()


library(Seurat)
library(ggplot2)
library(irGSEA)
seu_ec <- readRDS("~/seu_ec_715.rds")

seu_ec <- subset(seu_ec, Race %in% c("African", "Asian", "Caucasian"))
Idents(seu_ec) <- "Race"

seu_ec <- irGSEA.score(object = seu_ec, assay = "SCT", 
                       slot = "data", seeds = 123, ncores = 10,
                       custom = F, geneset = NULL, msigdb = T, 
                       species = "Homo sapiens", category = "H",  
                       subcategory = NULL, geneid = "symbol",
                       method = c("AUCell", "UCell", "singscore", 
                                  "ssgsea", "JASMINE", "viper"),
                       kcdf = 'Gaussian')
saveRDS(seu_ec, file = "seu_ec_race_score_H.rds")
result.dge.H <- irGSEA.integrate(object = seu_ec, 
                                 group.by = "Race",
                                 metadata = NULL, col.name = NULL,
                                 method = c("AUCell","UCell","singscore",
                                            "ssgsea", "JASMINE", "viper"))
saveRDS(result.dge.H, file = "result.dge.race.H.rds")

####
seu_ec <- irGSEA.score(object = seu_ec, assay = "SCT", 
                       slot = "data", seeds = 123, ncores = 10,
                       custom = F, geneset = NULL, msigdb = T, 
                       species = "Homo sapiens", category = "C2",  
                       subcategory = "CP:KEGG", geneid = "symbol",
                       method = c("AUCell", "UCell", "singscore", 
                                  "ssgsea", "JASMINE", "viper"),
                       kcdf = 'Gaussian')
saveRDS(seu_ec, file = "seu_ec_score_race.KEGG.rds")
result.dge.KEGG <- irGSEA.integrate(object = seu_ec, 
                                    group.by = "Race",
                                    metadata = NULL, col.name = NULL,
                                    method = c("AUCell","UCell","singscore",
                                               "ssgsea", "JASMINE", "viper"))
saveRDS(result.dge.KEGG, file = "result.dge.race.KEGG.rds")

pdf(file = "heatmap-race-H-RRA.pdf", width = 15, height = 15)
irGSEA.heatmap(object = result.dge.race.H, method = "RRA", top = 50, show.geneset = NULL, heatmap.width = 20, heatmap.heigh = 18)
dev.off()


pdf(file = "heatmap-race-KEGG-RRA.pdf", width = 15, height = 15)
irGSEA.heatmap(object = result.dge.race.KEGG, method = "RRA", top = 50, show.geneset = NULL, heatmap.width = 20, heatmap.heigh = 18)
dev.off()

pdf(file = "race-H-upset-RRA.pdf", width = 7, height = 5)
irGSEA.upset(object = result.dge.race.H, method = "RRA")
dev.off()


pdf(file = "race-kegg-upset-RRA.pdf", width = 7, height = 5)
irGSEA.upset(object = result.dge.race.KEGG, method = "RRA")
dev.off()


pdf(file = "race-H-upset-RRA.pdf", width = 7, height = 5)
irGSEA.upset(object = result.dge.race.H, method = "RRA")
dev.off()


pdf(file = "race-kegg-upset-RRA.pdf", width = 7, height = 5)
irGSEA.upset(object = result.dge.race.KEGG, method = "RRA")
dev.off()


pdf(file = "race-kegg-barplot-H.pdf", height = 7, width = 13)
irGSEA.barplot(object = result.dge.race.H, method = c("AUCell", "UCell", "singscore", "ssgsea", "JASMINE", "viper", "RRA"))
dev.off()


pdf(file = "race-kegg-barplot-KEGG.pdf", height = 7, width = 13)
irGSEA.barplot(object = result.dge.race.KEGG, method = c("AUCell", "UCell", "singscore", "ssgsea", "JASMINE", "viper", "RRA"))
dev.off()






