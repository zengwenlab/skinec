##Direder2022_GSE181316[seu1]
Direder2022_5 <- readRDS("~/Skin_project/Rds/Direder2022/Direder2022_5.rds")
table(Direder2022_5$orig.ident)
seu1 <- subset(Direder2022_5, orig.ident == "skin1")
table(seu1$orig.ident)

seu1$percent.mt <- NULL
seu1$doubFind_res <- NULL
seu1$doubFind_score <- NULL

rm(Direder2022_5);gc()

##Tabib2018[seu2]
Tabib2018_5 <- readRDS("~/Skin_project/Rds/Tabib2018/Tabib2018_5.rds")
table(Tabib2018_5$orig.ident)
seu2 <- Tabib2018_5
seu2$percent.mt <- NULL
seu2$doubFind_res <- NULL
seu2$doubFind_score <- NULL

rm(Tabib2018_5);gc()

##Vorstandlechner2021_GSE156326[seu3]
Vorstandlechner2021_5 <- readRDS("~/Skin_project/Rds/Vorstandlechner2021/\
Vorstandlechner2021_5.rds")
table(Vorstandlechner2021_5$orig.ident)
seu3 <- subset(Vorstandlechner2021_5, orig.ident %in% c("skin8", "skin9", "skin10"))
table(seu3$orig.ident)
seu3$orig.ident <- factor(seu3$orig.ident, levels = c("skin8", "skin9", "skin10"))
seu3$percent.mt <- NULL
seu3$doubFind_res <- NULL
seu3$doubFind_score <- NULL

rm(Vorstandlechner2021_5);gc()

###SoleBoldo2020_GSE130973[seu4]
SoleBoldo2020_5 <- readRDS("~/Skin_project/Rds/SoleBoldo2020/SoleBoldo2020_5.rds.rds")
table(SoleBoldo2020_5$orig.ident)
seu4 <- SoleBoldo2020_5
seu4$percent.mt <- NULL
seu4$doubFind_res <- NULL
seu4$doubFind_score <- NULL

rm(SoleBoldo2020_5);gc()

###Rustagi2022_GSE182208[seu5]
Rustagi2022_5 <- readRDS("~/Skin_project/Rds/Rustagi2022/Rustagi2022_5.rds")
table(Rustagi2022_5$orig.ident)
seu5 <- Rustagi2022_5
seu5$percent.mt <- NULL
seu5$doubFind_res <- NULL
seu5$doubFind_score <- NULL

rm(Rustagi2022_5);gc()

###Gao2021_GSE162183[seu6]
Gao2021_5 <- readRDS("~/Skin_project/Rds/Gao2021/Gao2021_5.rds")
table(Gao2021_5$orig.ident)
seu6 <- Gao2021_5
seu6$percent.mt <- NULL
seu6$doubFind_res <- NULL
seu6$doubFind_score <- NULL

rm(Gao2021_5);gc()

###Wang2021_GSE158924[seu7]
Wang2021_5 <- readRDS("~/Skin_project/Rds/Wang2021/Wang2021_5.rds")
table(Wang2021_5$orig.ident)
seu7 <- Wang2021_5
seu7$percent.mt <- NULL
seu7$doubFind_res <- NULL
seu7$doubFind_score <- NULL

rm(Wang2021_5);gc()

###Zou2021[seu8]
seu8 <- readRDS("~/Skin_project/Rds/Zou2021/Zou2021_5.rds")
table(seu8$orig.ident)
seu8$percent.mt <- NULL
seu8$doubFind_res <- NULL
seu8$doubFind_score <- NULL

###Ahlers2021[seu9]
seu9 <- readRDS("~/Skin_project/Rds/Ahlers2021/Ahlers2021_4.rds")
table(seu9$orig.ident)
seu9$percent.mt <- NULL
seu9$doubFind_res <- NULL
seu9$doubFind_score <- NULL

###Liu2022[seu10]
Liu2022_5 <- readRDS("~/Skin_project/Rds/Liu2022/Liu2022_5.rds")
table(Liu2022_5$orig.ident)
seu10 <- subset(Liu2022_5, orig.ident %in% c("skin43", "skin44", "skin45", "skin46"))
table(seu10$orig.ident)
seu10$percent.mt <- NULL
seu10$doubFind_res <- NULL
seu10$doubFind_score <- NULL

rm(Liu2022_5);gc()

###Xue2022[seu11]
seu11 <- readRDS("~/Skin_project/Rds/Xue2022/Xue2022_5.rds")
table(seu11$orig.ident)
seu11$percent.mt <- NULL
seu11$doubFind_res <- NULL
seu11$doubFind_score <- NULL

###GSE212447[seu12]
seu12 <- readRDS("~/Skin_project/Rds/GSE212447/GSE212447_5.rds")
seu12$percent.mt <- NULL
seu12$doubFind_res <- NULL
seu12$doubFind_score <- NULL

table(seu12$orig.ident)

###GSE202352[seu13]
seu13 <- readRDS("~/Skin_project/Rds/GSE202352/GSE202352_5.rds")
table(seu13$orig.ident)

rm(GSE202352_1);gc()

###GSE179162[seu14]
seu14 <- readRDS("~/Skin_project/Rds/GSE179162/GSE179162_5.rds")
table(seu14$orig.ident)


##########
seu <- merge(seu1, c(seu2, seu3, seu4, seu5, seu6, seu7, 
                     seu8, seu9, seu10, seu11, seu12, seu13, seu14))

seu$percent.mt <- NULL
seu$doubFind_res <- NULL
seu$doubFind_score <- NULL
table(seu$orig.ident)

seu$orig.ident <- factor(seu$orig.ident, 
                         levels = c("skin1", "skin2", "skin3", "skin4", "skin5", 
                                    "skin6", "skin7", "skin8", "skin9", "skin10", 
                                    "skin11", "skin12", "skin13", "skin14", "skin15", 
                                    "skin16", "skin17", "skin18", "skin19", "skin20", 
                                    "skin21", "skin22", "skin23", "skin24","skin25", 
                                    "skin26", "skin27", "skin28", "skin29", "skin30", 
                                    "skin31", "skin32", "skin33", "skin34", "skin35", 
                                    "skin36","skin37", "skin38", "skin39", "skin40", 
                                    "skin41", "skin42", "skin43", "skin44", "skin45", 
                                    "skin46", "skin47", "skin48","skin49", "skin50", 
                                    "skin51", "skin52", "skin53", "skin54", "skin55", 
                                    "skin56", "skin57", "skin58", "skin59", "skin60",
                                    "skin61", "skin62", "skin63", "skin64", "skin65"))

####
rm(counts_skin641, counts_skin642, counts_skin651, 
   counts_skin652, seu_list, seu1, seu2, seu3, seu4, 
   seu5, seu6, seu7, seu8, seu9, seu10, seu11, seu12, 
   seu13, seu14, seu641, seu642, seu651, seu652, skin64, skin65);gc()

rm(DoubletFinder_multisamples)

##########

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims = 1:30, reduction = "pca")
seu <- FindNeighbors(seu, dims = 1:30, reduction = "pca")
seu <- FindClusters(seu, resolution = 0.3)
DimPlot(seu, label = T, reduction = "umap", split.by = "orig.ident")


seu$Sample <- seu$orig.ident
table(seu$Sample)
data <- read_csv("info_table.csv")
data2 <- seu@meta.data
data2 <- rownames_to_column(data2,var = "barcodes")
metadata <- merge(data2,
                  data,
                  by.x='orig.ident',
                  by.y='Sample',
                  all=T)
metadata <- column_to_rownames(metadata,var = "barcodes")
seu <- AddMetaData(seu,metadata)



s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes)
seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")
seu <- SCTransform(seu, vars.to.regress = c("nCount_RNA", "percent.mt", 
                                            "S.Score", "G2M.Score"), verbose = T)

DefaultAssay(seu)
seu <- RunPCA(seu, verbose = T)

library(harmony)
seu <- RunHarmony(seu, group.by.vars = "orig.ident", max.iter.harmony = 20)

DefaultAssay(seu)
conflicts_prefer(SeuratObject::saveRDS)
saveRDS(seu, file = "seu_skin_SCT_harmony.rds")






