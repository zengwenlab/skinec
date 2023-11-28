seu_ec <- readRDS("~/skin_project/figure2ec/seu_ec_715.rds")

Idents(seu_ec)

seu_ec$Gender <- factor(seu_ec$Gender, levels = c("female", "male"))

seu_ec_gender <- subset(seu_ec, Gender %in% c("female", "male"))

pdf(file = "gender-umap.pdf", width = 13, height = 6)
DimPlot(seu_ec_gender, reduction = "umap", split.by = "Gender", label = T, 
        repel = T, group.by = "seurat_clusters")
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
pdf(file = "gender-propotion.pdf", width = 5, height = 10)
proportion_barplot(seu_ec_gender, cluster1 = "seurat_clusters", 
                   cluster2 = "Gender", col.cluster1 = cols_clusters)
dev.off()


seu_ec_gender$Sample <- droplevels(seu_ec_gender$Sample)

Cellratio <- prop.table(table(Idents(seu_ec_gender), seu_ec_gender$Sample), margin = 2)
Cellratio <- data.frame(Cellratio)
library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]

###添加分组信息
samples <- read_csv("1.csv")
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
  
  ###组间t检验分析
  labely = max(cellper_$percent)
  compare_means(percent ~ group$group,  data = cellper_)
  my_comparisons <- list( c("female", "male") )
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons, method = "t.test")
  pplist[[group_]] = pp1
}


pdf(file = "pro.pdf", width = 10, height = 15)
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

VolcanoPlot=function(dif, log2FC=log2(1.5), padj=0.05, 
                     label.symbols=NULL, label.max=30,
                     cols=c("#497aa2", "#ae3137"), title=""){
  if( all( !c("log2FoldChange", "padj", "symbol") %in% colnames(dif) )){
    stop("Colnames must include: log2FoldChange, padj, symbol")
  }
  rownames(dif)=dif$symbol
  
  # (1) define up and down
  dif$threshold="ns";
  dif[which(dif$log2FoldChange > log2FC & dif$padj <padj),]$threshold="up";
  dif[which(dif$log2FoldChange < (-log2FC) & dif$padj < padj),]$threshold="down";
  dif$threshold=factor(dif$threshold, levels=c('down','ns','up'))
  #head(dif)
  #
  tb2=table(dif$threshold); print(tb2)
  library(ggplot2)
  # (2) plot
  g1 = ggplot(data=dif, aes(x=log2FoldChange, y=-log10(padj), color=threshold)) +
    geom_point(alpha=0.8, size=0.8) +
    geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+
    geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
    labs(title= ifelse(""==title, "", paste("DEG:", title)))+
    xlab(bquote(Log[2]*FoldChange))+
    ylab(bquote(-Log[10]*italic(P.adj)) )+
    theme_classic(base_size = 14) +
    theme(legend.box = "horizontal",
          legend.position="top",
          legend.spacing.x = unit(0, 'pt'),
          legend.text = element_text( margin = margin(r = 20) ),
          legend.margin=margin(b= -10, unit = "pt"),
          plot.title = element_text(hjust = 0.5, size=10)
    ) +
    scale_color_manual('',labels=c(paste0("down(",tb2[[1]],')'),'ns',
                                   paste0("up(",tb2[[3]],')' )),
                       values=c(cols[1], "grey", cols[2]) )+
    guides(color=guide_legend(override.aes = list(size=3, alpha=1))); g1;
  # (3)label genes
  if(is.null(label.symbols)){
    dif.sig=dif[which(dif$threshold != "ns" ), ]
    len=nrow(dif.sig)
    if(len<label.max){
      label.symbols=rownames(dif.sig)
    }else{
      dif.sig=dif.sig[order(dif.sig$log2FoldChange), ]
      dif.sig= rbind(dif.sig[1:(label.max/2),], dif.sig[(len-label.max/2):len,])
      label.symbols=rownames(dif.sig)
    }
  }
  dd_text = dif[label.symbols, ]
  print((dd_text))
  # add text
  library(ggrepel)
  g1 + geom_text_repel(data=dd_text,
                       aes(x=log2FoldChange, y=-log10(padj), label=row.names(dd_text)),
                       #size=2.5, 
                       colour="black",alpha=1)
}


Idents(seu_ec_gender) <- "Gender"
deg_gender=FindMarkers(seu_ec_gender, ident.1 = "female", ident.2 = "male")


dif=data.frame(
  symbol=rownames(deg_gender),
  log2FoldChange=deg_gender$avg_log2FC,
  padj=deg_gender$p_val_adj
)

pdf(file = "femalevsmale.pdf", width = 6, height = 5)
VolcanoPlot(dif, padj=0.05, title="female vs male", 
            label.max = 10, cols=c("blue", "red"))
dev.off()


cluster5 <- subset(seu_ec_gender, seurat_clusters == "5")

Idents(cluster5) <- "Gender"
deg_gender5=FindMarkers(cluster5, ident.1 = "female", ident.2 = "male")

dif2=data.frame(
  symbol=rownames(deg_gender5),
  log2FoldChange=deg_gender5$avg_log2FC,
  padj=deg_gender5$p_val_adj
)

pdf(file = "femalevsmale-cluster5.pdf", width = 6, height = 5)
VolcanoPlot(dif2, padj=0.05, title="cluster5 female vs male", 
            label.max = 10, cols=c("blue", "red"))
dev.off()


cluster9 <- subset(seu_ec_gender, seurat_clusters == "9")

Idents(cluster9) <- "Gender"
deg_gender9=FindMarkers(cluster9, ident.1 = "female", ident.2 = "male")

dif2=data.frame(
  symbol=rownames(deg_gender9),
  log2FoldChange=deg_gender9$avg_log2FC,
  padj=deg_gender9$p_val_adj
)

pdf(file = "femalevsmale-cluster9.pdf", width = 6, height = 5)
VolcanoPlot(dif2, padj=0.05, title="cluster9 female vs male", 
            label.max = 10, cols=c("blue", "red"))
dev.off()

deg_gender$state <- ifelse(deg_gender$avg_log2FC>0, "up", "down")
deg_gender_up <- subset(deg_gender, state == "up")
deg_gender_down <- subset(deg_gender, state == "down")

#up
deg_gender_up <- rownames_to_column(deg_gender_up, var = "gene")
deg_gender_down <- rownames_to_column(deg_gender_down, var = "gene")


library(org.Hs.eg.db)
library(clusterProfiler)
go_subset_id <- bitr(deg_gender_up[,1],
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

pdf(file = "GO-all.pdf", width = 6, height = 8)
dotplot(GO)
dev.off()

library(org.Hs.eg.db)
library(clusterProfiler)
go_subset_id <- bitr(deg_gender_down[,1],
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

pdf(file = "GO-all-down.pdf", width = 6, height = 8)
dotplot(GO)
dev.off()


library(scMetabolism)
library(ggplot2)
library(rsvd)


seu_ec_gender <- sc.metabolism.Seurat(obj = seu_ec_gender, 
                                      method = "VISION", imputation = F, 
                                      ncores = 2, metabolism.type = "KEGG")

metabolism.matrix <- seu_ec_gender@assays$METABOLISM$score


DimPlot.metabolism(obj = seu_ec_gender, 
                   pathway = "Glycolysis/Gluconeogenesis", 
                   dimention.reduction.type = "umap", 
                   dimention.reduction.run = F, size = 1)

input.pathway <- rownames(seu_ec_gender@assays[["METABOLISM"]][["score"]])[1:30]

seu_ec_gender$seurat_clusters <- factor(seu_ec_gender$seurat_clusters, 
                                        levels = c("0", "1", "2", "3", "4", 
                                                   "5", "6", "7", "8", "9", "10"))

pdf(file = "cluster_metabolism.pdf", height = 10, width = 10)
DotPlot.metabolism(obj = seu_ec_gender, 
                   pathway = input.pathway, 
                   phenotype = "seurat_clusters", 
                   norm = "y")
dev.off()

pdf(file = "gender_metabolism.pdf", height = 10, width = 5)
DotPlot.metabolism(obj = seu_ec_gender, 
                   pathway = input.pathway, 
                   phenotype = "Gender", 
                   norm = "y")
dev.off()


sce_Metal_exp = seu_ec_gender
sce_Metal_exp$celltype = sce_Metal_exp$seurat_clusters
mscore_data = data.frame(t(sce_Metal_exp@assays[["METABOLISM"]][["score"]]),sce_Metal_exp$celltype)
avg_sM=aggregate(mscore_data[,1:ncol(mscore_data)-1],list(mscore_data$sce_Metal_exp.celltype),mean)
rownames(avg_sM) = avg_sM$Group.1
avg_sM=data.frame(t(avg_sM[,-1]))
avg_sM$KEGG = rownames(sce_Metal_exp@assays[["METABOLISM"]][["score"]])
rownames(avg_sM)=avg_sM$KEGG

c_k_l = c()
for(c in c(1:ncol(avg_sM))){
  c_k=avg_sM[order(avg_sM[,c]),]$KEGG[1:5]
  c_k_l=c(c_k_l,c_k)
}
c_k_l= unique(c_k_l)
c_k_d = avg_sM[avg_sM$KEGG %in%c_k_l,]


rownames(c_k_d) = c_k_d$KEGG

pdf(file = "kegg_metabolism.pdf", width = 10, height = 10)
pheatmap::pheatmap(c_k_d[,-ncol(c_k_d)],show_colnames = T,scale='row')
dev.off()


countexp.Seurat <- seu_ec_gender
input.pathway <- rownames(seu_ec_gender@assays[["METABOLISM"]][["score"]])[18:20]
pdf(file = "boxplot.metabolism2.pdf", width = 12, height = 4)
BoxPlot.metabolism(obj = countexp.Seurat, 
                   pathway = input.pathway, 
                   phenotype = "seurat_clusters", ncol = 3)
dev.off()


countexp.Seurat <- seu_ec_gender
input.pathway <- rownames(seu_ec_gender@assays[["METABOLISM"]][["score"]])[18:20]
pdf(file = "boxplot.metabolism.pdf", width = 12, height = 4)
BoxPlot.metabolism(obj = countexp.Seurat, 
                   pathway = input.pathway, 
                   phenotype = "Gender", ncol = 3)
dev.off()


sce_Metal_exp = seu_ec_gender
sce_Metal_exp$celltype = sce_Metal_exp$Gender
mscore_data = data.frame(t(sce_Metal_exp@assays[["METABOLISM"]][["score"]]),sce_Metal_exp$celltype)
avg_sM=aggregate(mscore_data[,1:ncol(mscore_data)-1],list(mscore_data$sce_Metal_exp.celltype),mean)
rownames(avg_sM) = avg_sM$Group.1
avg_sM=data.frame(t(avg_sM[,-1]))
avg_sM$KEGG = rownames(sce_Metal_exp@assays[["METABOLISM"]][["score"]])
rownames(avg_sM)=avg_sM$KEGG

c_k_l = c()
for(c in c(1:ncol(avg_sM))){
  c_k=avg_sM[order(avg_sM[,c]),]$KEGG[1:10]
  c_k_l=c(c_k_l,c_k)
}
c_k_l= unique(c_k_l)
c_k_d = avg_sM[avg_sM$KEGG %in%c_k_l,]


rownames(c_k_d) = c_k_d$KEGG

pdf(file = "kegg_metabolism2.pdf", width = 5, height = 8)
pheatmap::pheatmap(c_k_d[,-ncol(c_k_d)],show_colnames = T,scale='row')
dev.off()


cluster5 <- sc.metabolism.Seurat(obj = cluster5, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG")
input.pathway <- rownames(cluster5@assays[["METABOLISM"]][["score"]])[1:50]
pdf(file = "gender_cluster5_metabolism.pdf", height = 10, width = 5)
DotPlot.metabolism(obj = cluster5, 
                   pathway = input.pathway, 
                   phenotype = "Gender", 
                   norm = "y")
dev.off()


sce_Metal_exp = cluster5
sce_Metal_exp$celltype = sce_Metal_exp$Gender
mscore_data = data.frame(t(sce_Metal_exp@assays[["METABOLISM"]][["score"]]),sce_Metal_exp$celltype)
avg_sM=aggregate(mscore_data[,1:ncol(mscore_data)-1],list(mscore_data$sce_Metal_exp.celltype),mean)
rownames(avg_sM) = avg_sM$Group.1
avg_sM=data.frame(t(avg_sM[,-1]))
avg_sM$KEGG = rownames(sce_Metal_exp@assays[["METABOLISM"]][["score"]])
rownames(avg_sM)=avg_sM$KEGG

c_k_l = c()
for(c in c(1:ncol(avg_sM))){
  c_k=avg_sM[order(avg_sM[,c]),]$KEGG[1:10]
  c_k_l=c(c_k_l,c_k)
}
c_k_l= unique(c_k_l)
c_k_d = avg_sM[avg_sM$KEGG %in%c_k_l,]


rownames(c_k_d) = c_k_d$KEGG

pdf(file = "kegg_metabolism_cluster5.pdf", width = 8, height = 8)
pheatmap::pheatmap(c_k_d[,-ncol(c_k_d)],show_colnames = T,scale='row')
dev.off()


group1 <- subset(seu_ec_gender, Sample %in% c("skin9", "skin1"))


table(group1$orig.ident)
group1$orig.ident <- droplevels(group1$orig.ident)
group1$Sample <- droplevels(group1$Sample)
table(group1$Gender)

Idents(group1)
degs_group1 <- FindMarkers(group1, ident.1 = "female", ident.2 = "male")

degs_group1 <- rownames_to_column(degs_group1, var = "gene")
write_csv(degs_group1, file = "degs_group1.csv")

group2 <- subset(seu_ec_gender, Sample %in% c("skin49", "skin47"))


table(group2$orig.ident)
group2$orig.ident <- droplevels(group2$orig.ident)
group2$Sample <- droplevels(group2$Sample)
table(group2$Gender)

Idents(group2)
degs_group2 <- FindMarkers(group2, ident.1 = "female", ident.2 = "male")

degs_group2 <- rownames_to_column(degs_group2, var = "gene")
write_csv(degs_group2, file = "degs_group2.csv")

####skin38 vs. skin18

group22 <- subset(seu_ec_gender, Sample %in% c("skin38", "skin18"))
table(group22$orig.ident)
group22$orig.ident <- droplevels(group22$orig.ident)
group22$Sample <- droplevels(group22$Sample)
table(group22$Gender)

Idents(group22)
degs_group22 <- FindMarkers(group22, ident.1 = "female", ident.2 = "male")

degs_group22 <- rownames_to_column(degs_group22, var = "gene")
write_csv(degs_group22, file = "degs_group22.csv")


group3 <- subset(seu_ec_gender, Sample %in% c("skin53", "skin51"))

table(group3$orig.ident)
group3$orig.ident <- droplevels(group3$orig.ident)
group3$Sample <- droplevels(group3$Sample)
table(group3$Gender)

Idents(group3)
degs_group3 <- FindMarkers(group3, ident.1 = "female", ident.2 = "male")

degs_group3 <- rownames_to_column(degs_group3, var = "gene")
write_csv(degs_group3, file = "degs_group3.csv")

degs_group1$label <- ifelse(degs_group1$avg_log2FC>0, "up", "down")
degs_group2$label <- ifelse(degs_group2$avg_log2FC>0, "up", "down")
degs_group3$label <- ifelse(degs_group3$avg_log2FC>0, "up", "down")

degs_group1_up <- subset(degs_group1, label=="up")
degs_group2_up <- subset(degs_group2, label=="up")
degs_group3_up <- subset(degs_group3, label=="up")

degs_group1_down <- subset(degs_group1, label=="down")
degs_group2_down <- subset(degs_group2, label=="down")
degs_group3_down <- subset(degs_group3, label=="down")


library(venn)
library(VennDiagram)

pdf(file = "venn-up.pdf", width = 8, height = 8)
venn(list(group1=degs_group1_up$gene, group2=degs_group2_up$gene, group3=degs_group3_up$gene), zcolor = 2:4, box=F, sncs = 2, ilcs = 1.5)
dev.off()


pdf(file = "venn-down.pdf", width = 8, height = 8)
venn(list(group1=degs_group1_down$gene, group2=degs_group2_down$gene, group3=degs_group3_down$gene), zcolor = 2:4, box=F, sncs = 2, ilcs = 1.5)
dev.off()


df_inter <- get.venn.partitions(list(group1=degs_group1_up$gene, group2=degs_group2_up$gene, group3=degs_group3_up$gene))
for(i in 1:nrow(df_inter)){
  df_inter[i, "values"] <- paste(df_inter[[i, "..values.."]], collapse = ",")
}
df_inter <- df_inter[-c(4,5)]
write_csv(df_inter, file = "venn-up.csv")


df_inter <- get.venn.partitions(list(group1=degs_group1_down$gene, group2=degs_group2_down$gene, group3=degs_group3_down$gene))
for(i in 1:nrow(df_inter)){
  df_inter[i, "values"] <- paste(df_inter[[i, "..values.."]], collapse = ",")
}
df_inter <- df_inter[-c(4,5)]
write_csv(df_inter, file = "venn-down.csv")


upgene <- read_csv(file = "venn-up.csv")

upgene <- upgene[1,5]
upgene <- str_split(upgene, ",")
upgene <- as.data.frame(upgene)
colnames(upgene)[1] <- "SYMBOL"

library(clusterProfiler)
library(org.Hs.eg.db)
go_subset_id <- bitr(upgene[,1],
                     fromType="SYMBOL",
                     toType=c("ENTREZID","ENSEMBL","SYMBOL"),
                     OrgDb="org.Hs.eg.db")

GO_Result <- enrichGO(go_subset_id$ENTREZID,
                      OrgDb = org.Hs.eg.db, 
                      ont='ALL',
                      pAdjustMethod = 'BH',
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.2,
                      keyType = 'ENTREZID')

pdf(file = "go-up.pdf", width = 7, height = 5)
barplot(GO_Result)  
dev.off()



downgene <- read_csv(file = "venn-down.csv")


downgene <- downgene[1,5]
downgene <- str_split(downgene, ",")
downgene <- as.data.frame(downgene)
colnames(downgene)[1] <- "SYMBOL"

library(clusterProfiler)
library(org.Hs.eg.db)
go_subset_id <- bitr(downgene[,1],
                     fromType="SYMBOL",
                     toType=c("ENTREZID","ENSEMBL","SYMBOL"),
                     OrgDb="org.Hs.eg.db")

GO_Result <- enrichGO(go_subset_id$ENTREZID,
                      OrgDb = org.Hs.eg.db, 
                      ont='ALL',
                      pAdjustMethod = 'BH',
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.2,
                      keyType = 'ENTREZID')

pdf(file = "go-downgene.pdf", width = 7, height = 5)
barplot(GO_Result)  
dev.off()



group_seu_ec <- subset(seu_ec_gender, Sample %in% c("skin9", "skin1", "skin49", "skin47", "skin53", "skin51"))
group_seu_ec$orig.ident <- droplevels(group_seu_ec$orig.ident)
group_seu_ec$Sample <- droplevels(group_seu_ec$Sample)
Idents(group_seu_ec)
VlnPlot(group_seu_ec, features = upgene$SYMBOL, pt.size = 0)
pdf(file = "upgene.pdf", width = 10, height = 10)
VlnPlot(group_seu_ec, features = c("C6orf48", "RPL31", "EVA1C", "RPL23A", "SOCS2", "RND1", "IFIT3", "KLF7", "BAG3"), ncol = 3, pt.size = 0)
dev.off()


pdf(file = "downgene.pdf", width = 10, height = 10)
VlnPlot(group_seu_ec, features = c("CAVIN1", "SNHG7", "CD81", "ELOB", "PLCG2", "KRTCAP2", "SET", "ETS1", "CAST"), ncol = 3, pt.size = 0)
dev.off()

library(Seurat)
library(ggplot2)
library(irGSEA)
seu_ec <- readRDS("~/Skin_project/seu_ec_715.rds")

seu_ec <- subset(seu_ec, Gender %in% c("female", "male"))
seu_ec$Gender <- factor(seu_ec$Gender, levels = c("female", "male"))


Idents(seu_ec) <- "Gender"

###Hallmark
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

####KEGG
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


pdf(file = "heatmap-H-RRA.pdf", width = 9, height = 6)
irGSEA.heatmap(object = result.dge.gender.H, method = "RRA", top = 50, show.geneset = NULL, heatmap.width = 20, heatmap.heigh = 18)
dev.off()

pdf(file = "heatmap-KEGG-RRA.pdf", width = 9, height = 6)
irGSEA.heatmap(object = result.dge.gender.KEGG, method = "RRA", top = 50, show.geneset = NULL, heatmap.width = 20, heatmap.heigh = 18)
dev.off()










