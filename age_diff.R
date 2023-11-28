seu_ec_715 <- readRDS("~/skin_project/figure2ec/seu_ec_715.rds")
library(Seurat)
library(ggplot2)

seu_ec <- seu_ec_715
rm(seu_ec_715);gc()

DimPlot(seu_ec, group.by = "seurat_clusters", split.by = "Age2")


seu_ec_age2 <- subset(seu_ec, Age2 %in% c("Y", "M", "O"))

pdf(file = "umap-age2.pdf", height = 5, width = 10)
DimPlot(seu_ec_age2, group.by = "seurat_clusters", split.by = "Age2", label = T, repel = T)
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
pdf(file = "Age2-propotion.pdf", width = 5, height = 8)
proportion_barplot(seu_ec_age2, cluster1 = "seurat_clusters", 
                   cluster2 = "Age2", col.cluster1 = cols_clusters)
dev.off()

seu_ec_age2 <- subset(seu_ec, Age2 %in% c("Y", "M", "O"))
rm(seu_ec);gc()
seu_ec_age2$orig.ident <- droplevels(seu_ec_age2$orig.ident)
seu_ec_age2$Sample <- droplevels(seu_ec_age2$Sample)
seu_ec_age2$Age2 <- droplevels(seu_ec_age2$Age2)



Cellratio <- prop.table(table(Idents(seu_ec_age2), seu_ec_age2$Sample), margin = 2)
Cellratio <- data.frame(Cellratio)
library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]

###添加分组信息
samples <- read_csv("11.csv")
colnames(samples) <- c("sample", "group")
rownames(samples) <- samples$sample

cellper$sample <- samples[rownames(cellper),'sample']#R添加列
cellper$group <- samples[rownames(cellper),'group']#R添加列
cellper <- as.data.frame(cellper)


###作图展示
pplist = list()
sce_groups = c('0','1', '2', '3', '4', '5', '6', '7', '8', '9', '10')
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflicts_prefer(dplyr::select)


for(group_ in sce_groups){
  cellper_  = cellper %>% select(one_of(c('sample','group', group_)))#选择一组数据
  colnames(cellper_) = c('sample','group','percent')#对选择数据列命名
  cellper_$percent = as.numeric(cellper_$percent)#数值型数据
  cellper_ <- cellper_ %>% group_by(group$group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                            lower = quantile(percent, 0.25),
                                                            mean = mean(percent),
                                                            median = median(percent))#上下分位数
  print(group_)
  print(cellper_$median)
  
  cellper_$`group$group` <- factor(cellper_$`group$group`, levels = c("Y", "M", "O"))
  
  pp1 = ggplot(cellper_, aes(x=group$group, y=percent)) + #ggplot作图
    geom_jitter(shape = 21,aes(fill=group$group),width = 0.2, size=5) + 
    stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          plot.title = element_text(size = 10,face = 'plain'),
          legend.position = 'none') + 
    labs(title = group_,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  0.25)
  
  
  ###组间检验分析
  labely = max(cellper_$percent)
  compare_means(percent ~ group$group,  data = cellper_)
  my_comparisons <- list( c("Y", "M"), c("Y", "O"), c("M", "O") )
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons, method = "anova") 
  pplist[[group_]] = pp1
}


pdf(file = "pro22.pdf", width = 10, height = 12)
plot_grid(pplist[['0']],
          pplist[['1']],
          pplist[['2']],
          pplist[['3']],
          pplist[['4']],
          pplist[['5']],
          pplist[['6']],
          pplist[['7']],
          pplist[['8']],
          pplist[['9']],
          pplist[['10']], ncol = 3)
dev.off()


gene_age2 <- AverageExpression(seu_ec_age2, group.by = "Age2", assays = "SCT")
gene_age2 <- as.data.frame(gene_age2)
write.csv(gene_age2, file = "gene_age2.csv")


gene7 <- c('ADAMTS4', 'ADAMTS9', 'CXCL2', 'HLA-DRA', 'HLA-DRB1', 
           'GJA1', 'TNFAIP3', 'SOD2', 'IL6', 'AQP1', 'RPS20', 'NNMT', 
           'EMP1', 'CSF3', 'VMP1', 'DDX5', 'TCF4', 'ADAMTS1', 'MT-ND2', 'MT-CO1', 'MT-CO2', 'MT-CO3', 'MT-ND5')


library(org.Hs.eg.db)
library(clusterProfiler)
go_subset_id <- bitr(gene7,
                     fromType="SYMBOL",
                     toType=c("ENTREZID","ENSEMBL","SYMBOL"),
                     OrgDb="org.Hs.eg.db")
library(dplyr)
GO=enrichGO(gene = go_subset_id$ENTREZID,
            OrgDb = org.Hs.eg.db, 
            pvalueCutoff =1,	
            qvalueCutoff = 1,
            ont="all",	
            readable =T)

pdf(file = "GO-alldown-bar.pdf", width = 8, height = 6)
barplot(GO)
dev.off()

pdf(file = "gene7.pdf", height = 5, width = 8)
VlnPlot(seu_ec_age2, features = c('CXCL2', 'SOD2', 'IL6', 'CSF3','VMP1', 'NNMT'), ncol = 3, pt.size = 0)
dev.off()

Idents(seu_ec_age2) <- "Age2"

degs_age2 <- FindAllMarkers(seu_ec_age2, only.pos = T)

library(dplyr)
Degs_sig <- subset(degs_age2, p_val_adj < 0.01)
Degs_sig %>% group_by(gene) %>% top_n(1, avg_log2FC) -> Degs_sig_no_dup
Degs_sig_no_dup  %>% group_by(cluster) %>% top_n(5, -p_val_adj) -> top5
top5 %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

library(scRNAtoolVis)
pdf(file = "heatmap.ec-age2.pdf", width = 5, height = 5)
AverageHeatmap(object = seu_ec_age2, group.by = "Age2", annoCol = T, myanCol = cols_clusters[2:4],
               markerGene = Degs_sig_no_dup$gene, row_title = "", 
               showRowNames = F, markGenes = top5$gene)
dev.off()


library(monocle)
DefaultAssay(seu_ec_age2)
Idents(seu_ec_age2)

degs_age2 <- FindAllMarkers(seu_ec_age2)

data <- as(as.matrix(seu_ec_age2@assays$SCT@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seu_ec_age2@meta.data)
gene_annotation=data.frame(gene_short_name = rownames(seu_ec_age2[["SCT"]]),
                           stringsAsFactors=F)
rownames(gene_annotation)<-gene_annotation$gene_short_name
fd <- new("AnnotatedDataFrame", data = gene_annotation)

HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       expressionFamily=negbinomial.size())
library(dplyr)
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)


HSMM_ordering_genes=unique(degs_age2$gene)
HSMM <- setOrderingFilter(HSMM,ordering_genes = HSMM_ordering_genes) 

HSMM <- reduceDimension(HSMM,
                        max_components = 2,
                        reduction_method = "DDRTree") 
HSMM <- orderCells(HSMM)


cols_Set1 <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
plot_cell_trajectory(HSMM, color_by = "Age2", show_branch_points = T) + scale_color_manual(values = cols_Set1)


plot_cell_trajectory(HSMM, color_by = "Pseudotime")


HSMM2 <- reduceDimension(HSMM,
                         max_components = 3,
                         reduction_method = "DDRTree") 
HSMM2 <- orderCells(HSMM2)
plot_cell_trajectory(HSMM2, color_by = "Age2", show_branch_points = T) + scale_color_manual(values = cols_Set1)


diff_test_res <- differentialGeneTest(HSMM,fullModelFormulaStr = "~ Age2")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

cds_DGT <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(cds_DGT)
cds_DGT <- reduceDimension(cds_DGT, max_components = 2,reduction_method = 'DDRTree')
cds_DGT <- orderCells(cds_DGT)
#cds_DGT <- orderCells(cds_DGT, root_state = NULL, num_paths = NULL, reverse = T) 


plot_cell_trajectory(cds_DGT, color_by = "Age2", show_branch_points = T) + scale_color_manual(values = cols_Set1)


diff_test_res <- differentialGeneTest(HSMM,fullModelFormulaStr = "~Age2")
ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:2000]
HSMM <- setOrderingFilter(HSMM, ordering_genes)
HSMM <- reduceDimension(HSMM,
                        max_components = 2,
                        reduction_method = "DDRTree") 
HSMM <- orderCells(HSMM)
plot_cell_trajectory(HSMM, color_by = "Age2", show_branch_points = T) + scale_color_manual(values = cols_Set1)

plot_cell_trajectory(HSMM, color_by = "Pseudotime")
plot_cell_trajectory(HSMM, color_by = "State", show_branch_points = T) + scale_color_manual(values = cols_Set1)
HSMM <- orderCells(HSMM, root_state = 3, num_paths = 1, reverse = T)

pdf(file = "pseudotime.pdf", width = 8, height = 6)
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
dev.off()

pdf(file = "State.pdf", width = 8, height = 6)
plot_cell_trajectory(HSMM, color_by = "State", show_branch_points = T) + scale_color_manual(values = cols_Set1)
dev.off()

pdf(file = "Age2-gender.pdf", width = 10, height = 6)
plot_cell_trajectory(HSMM, color_by = "Age2", 
                     show_branch_points = T) + scale_color_manual(values = cols_Set1) + 
  facet_wrap(~Gender, nrow = 1)
dev.off()

cds_DGT_pseudotimegenes <- differentialGeneTest(HSMM, fullModelFormulaStr = "~sm.ns(Pseudotime)")
cds_DGT_pseudotimegenes_sig <- subset(cds_DGT_pseudotimegenes, qval < 0.01)
use_genes <- row.names(cds_DGT_pseudotimegenes)[order(cds_DGT_pseudotimegenes$qval)][1:2000]


saveRDS(cds_DGT_pseudotimegenes, file = "cds_DGT_pseudotimegenes.rds")


library(dplyr)
Time_genes <- use_genes %>% as.character()
p3I <- plot_pseudotime_heatmap(HSMM[Time_genes,], 
                               num_cluster = 4, 
                               show_rownames = T, 
                               return_heatmap = T)
ggsave(p3I, filename = "Figure3I.pdf", height = 20, width = 7)

p3I$tree_row

clusters <- cutree(p3I$tree_row, k = 4)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)


cluster1 <- subset(clustering, Gene_Clusters == "1")
cluster2 <- subset(clustering, Gene_Clusters == "2")
cluster3 <- subset(clustering, Gene_Clusters == "3")
cluster4 <- subset(clustering, Gene_Clusters == "4")

library(dplyr)
clustering <- tibble::rownames_to_column(clustering, "gene")
write.csv(clustering, file = "cluster.gene.csv")



library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(ggplot2)
library(stringr)
library(openxlsx)
library(org.Hs.eg.db)

######cluster1

entrezid_1 = mapIds(x = org.Hs.eg.db,
                    keys = rownames(cluster1),
                    keytype = "SYMBOL",
                    column = "ENTREZID")

entrezid_1  = na.omit(entrezid_1)
entrezid_1 = data.frame(entrezid_1)
head(entrezid_1)



GO_enrich = enrichGO(gene = entrezid_1[,1],
                     OrgDb = org.Hs.eg.db, 
                     keyType = "ENTREZID", 
                     ont = "BP", 
                     pvalueCutoff = 1,
                     qvalueCutoff = 1, 
                     readable = T) 
GO_enrich_1  = data.frame(GO_enrich) 
write.csv(GO_enrich_1,'GO_enrich_1.csv')

#####cluster2
entrezid_1 = mapIds(x = org.Hs.eg.db,
                    keys = rownames(cluster2),
                    keytype = "SYMBOL",
                    column = "ENTREZID")

entrezid_1  = na.omit(entrezid_1)
entrezid_1 = data.frame(entrezid_1)
head(entrezid_1)



GO_enrich = enrichGO(gene = entrezid_1[,1],
                     OrgDb = org.Hs.eg.db, 
                     keyType = "ENTREZID", 
                     ont = "BP", 
                     pvalueCutoff = 1,
                     qvalueCutoff = 1, 
                     readable = T) 
GO_enrich_1  = data.frame(GO_enrich) 
write.csv(GO_enrich_1,'GO_enrich_2.csv')



#####cluster3
entrezid_1 = mapIds(x = org.Hs.eg.db,
                    keys = rownames(cluster3),
                    keytype = "SYMBOL",
                    column = "ENTREZID")

entrezid_1  = na.omit(entrezid_1)
entrezid_1 = data.frame(entrezid_1)
head(entrezid_1)



GO_enrich = enrichGO(gene = entrezid_1[,1],
                     OrgDb = org.Hs.eg.db, 
                     keyType = "ENTREZID", 
                     ont = "BP", 
                     pvalueCutoff = 1,
                     qvalueCutoff = 1, 
                     readable = T) 
GO_enrich_1  = data.frame(GO_enrich) 
write.csv(GO_enrich_1,'GO_enrich_3.csv')

#####cluster4
entrezid_1 = mapIds(x = org.Hs.eg.db,
                    keys = rownames(cluster4),
                    keytype = "SYMBOL",
                    column = "ENTREZID")

entrezid_1  = na.omit(entrezid_1)
entrezid_1 = data.frame(entrezid_1)
head(entrezid_1)



GO_enrich = enrichGO(gene = entrezid_1[,1],
                     OrgDb = org.Hs.eg.db, 
                     keyType = "ENTREZID", 
                     ont = "BP", 
                     pvalueCutoff = 1,
                     qvalueCutoff = 1, 
                     readable = T) 
GO_enrich_1  = data.frame(GO_enrich) 
write.csv(GO_enrich_1,'GO_enrich_4.csv')


library(Seurat)
library(ggplot2)
library(irGSEA)
seu_ec <- readRDS("~/skin_project/figure2ec/seu_ec_715.rds")

seu_ec <- subset(seu_ec, Age2 %in% c("Y", "M", "O"))
Idents(seu_ec) <- "Age2"
seu_ec$Age2 <- droplevels(seu_ec$Age2)

seu_ec <- irGSEA.score(object = seu_ec, assay = "SCT", 
                       slot = "data", seeds = 123, ncores = 10,
                       custom = F, geneset = NULL, msigdb = T, 
                       species = "Homo sapiens", category = "H",  
                       subcategory = NULL, geneid = "symbol",
                       method = c("AUCell", "UCell", "singscore", 
                                  "ssgsea", "JASMINE", "viper"),
                       kcdf = 'Gaussian')
saveRDS(seu_ec, file = "seu_ec_gender_score_H.rds")
result.dge.H <- irGSEA.integrate(object = seu_ec, 
                                 group.by = "Gender",
                                 metadata = NULL, col.name = NULL,
                                 method = c("AUCell","UCell","singscore",
                                            "ssgsea", "JASMINE", "viper"))
saveRDS(result.dge.H, file = "result.dge.gender.H.rds")

####
seu_ec <- irGSEA.score(object = seu_ec, assay = "SCT", 
                       slot = "data", seeds = 123, ncores = 10,
                       custom = F, geneset = NULL, msigdb = T, 
                       species = "Homo sapiens", category = "C2",  
                       subcategory = "CP:KEGG", geneid = "symbol",
                       method = c("AUCell", "UCell", "singscore", 
                                  "ssgsea", "JASMINE", "viper"),
                       kcdf = 'Gaussian')
saveRDS(seu_ec, file = "seu_ec_score_gender.KEGG.rds")
result.dge.KEGG <- irGSEA.integrate(object = seu_ec, 
                                    group.by = "Gender",
                                    metadata = NULL, col.name = NULL,
                                    method = c("AUCell","UCell","singscore",
                                               "ssgsea", "JASMINE", "viper"))
saveRDS(result.dge.KEGG, file = "result.dge.gender.KEGG.rds")

pdf(file = "heatmap.age.H.pdf", width = 10, height = 10)
irGSEA.heatmap(object = result.dge.age.H, 
               method = "RRA",
               top = 50, 
               cluster.levels = c("Y", "M", "O"))

dev.off()



pdf(file = "heatmap.age.KEGG.pdf", width = 10, height = 10)
irGSEA.heatmap(object = result.dge.age.KEGG, 
               method = "RRA",
               top = 50,
               cluster.levels = c("Y", "M", "O"))

dev.off()
pdf(file = "upset.age.H.pdf", width = 10, height = 10)
irGSEA.upset(object = result.dge.age.H, 
             method = "RRA")
dev.off()



pdf(file = "upset.age.KEGG.pdf", width = 10, height = 10)
irGSEA.upset(object = result.dge.age.KEGG, 
             method = "RRA")
dev.off()




pdf(file = "barplot.age.H.pdf", width = 10, height = 10)
irGSEA.barplot(object = result.dge.age.H,
               method = c("AUCell", "UCell", "singscore",
                          "ssgsea", "JASMINE", "viper", "RRA"))
dev.off()


pdf(file = "barplot.age.KEGG.pdf", width = 10, height = 10)
irGSEA.barplot(object = result.dge.age.KEGG,
               method = c("AUCell", "UCell", "singscore",
                          "ssgsea", "JASMINE", "viper", "RRA"))
dev.off()

pdf(file = "ANGIOGENESIS.pdf", height = 10, width = 10)
irGSEA.halfvlnplot(object = seu_ec_age_score_H,
                   method = "AUCell",
                   show.geneset = "HALLMARK-ANGIOGENESIS")
dev.off()




