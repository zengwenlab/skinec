
library(CellChat)
library(Seurat)
table(seu$celltype)

seu$celltype1 <- ifelse((seu$celltype %in% c("Epi1","Epi2")), "Epi", 
                        ifelse(seu$celltype %in% c("Fib1","Fib2","Fib3"), "Fib", 
                               ifelse(seu$celltype %in% c("Smc1","Smc2"), "Smc",
                                      ifelse(seu$celltype %in% c("Imm1","Imm2"), "Imm", as.character(seu$celltype)))))
table(seu$celltype1)

seu$celltype1 <- factor(seu$celltype1, levels = c("Epi", "Fib", "Imm", "Smc", "VEc", "Melan", "LEc", "Swan", "Mac", "Gland", "Red", "Mast"))
Idents(seu) <- seu$celltype1
seu <- subset(seu, celltype1!= "Red")
seu$celltype1 <- droplevels(seu$celltype1)
Idents(seu) <- seu$celltype1

female <- subset(seu, Gender == "female")
DefaultAssay(female) <- "RNA"
female$S.Score <- NULL
female$SCT_snn_res.0.45 <- NULL
female$SCT_snn_res.0.4 <- NULL
female$SCT_snn_res.0.35 <- NULL
female$SCT_snn_res.0.3 <- NULL
female$SCT_snn_res.0.25 <- NULL
female$SCT_snn_res.0.2 <- NULL
female$SCT_snn_res.0.15 <- NULL
female$SCT_snn_res.0.1 <- NULL
female$SCT_snn_res.0.05 <- NULL
female$SCT_snn_res.0 <- NULL
female$nFeature_SCT <- NULL
female$percent.mt <- NULL
female$nCount_SCT <- NULL
female$Phase <- NULL
female$G2M.Score <- NULL
female$seurat_clusters <- NULL
female@assays$SCT <- NULL
female@reductions$pca <- NULL
female@reductions$harmony <- NULL
female@reductions$umap <- NULL


male <- subset(seu, Gender == "male")
DefaultAssay(male) <- "RNA"
male$S.Score <- NULL
male$SCT_snn_res.0.45 <- NULL
male$SCT_snn_res.0.4 <- NULL
male$SCT_snn_res.0.35 <- NULL
male$SCT_snn_res.0.3 <- NULL
male$SCT_snn_res.0.25 <- NULL
male$SCT_snn_res.0.2 <- NULL
male$SCT_snn_res.0.15 <- NULL
male$SCT_snn_res.0.1 <- NULL
male$SCT_snn_res.0.05 <- NULL
male$SCT_snn_res.0 <- NULL
male$nFeature_SCT <- NULL
male$percent.mt <- NULL
male$nCount_SCT <- NULL
male$Phase <- NULL
male$G2M.Score <- NULL
male$seurat_clusters <- NULL
male@assays$SCT <- NULL
male@reductions$pca <- NULL
male@reductions$harmony <- NULL
male@reductions$umap <- NULL

options(stringsAsFactors = FALSE)
options(futrue.globlas.Maxsize=2000*1024**3)
conflicts_prefer(SeuratObject::saveRDS)



###female
data.input.female <- GetAssayData(female, assay = "RNA", slot = "data")
labels <- Idents(female)
identity.female <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.female <- createCellChat(object = data.input.female)
cellchat.seu.female <- addMeta(cellchat.seu.female, meta = identity.female, 
                               meta.name = "labels")
cellchat.seu.female <- setIdent(cellchat.seu.female, ident.use = "labels") 

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.female@DB <- CellChatDB.use

cellchat.seu.female <- subsetData(cellchat.seu.female)
cellchat.seu.female <- identifyOverExpressedGenes(cellchat.seu.female)
cellchat.seu.female <- identifyOverExpressedInteractions(cellchat.seu.female)
cellchat.seu.female <- projectData(cellchat.seu.female, PPI.human)
cellchat.seu.female <- computeCommunProb(cellchat.seu.female)
cellchat.seu.female <- filterCommunication(cellchat.seu.female, min.cells = 10)
df.net.female <- subsetCommunication(cellchat.seu.female, slot.name = "net")
write.csv(df.net.female, "net.LR.female.csv")


cellchat.seu.female <- computeCommunProbPathway(cellchat.seu.female)
df.netp.female <- subsetCommunication(cellchat.seu.female, slot.name = "netP")
write.csv(df.netp.female, "net.patway.female.csv")

cellchat.seu.female <- aggregateNet(cellchat.seu.female)
cellchat.seu.female <- netAnalysis_computeCentrality(cellchat.seu.female)
saveRDS(cellchat.seu.female, file = "cellchat.female.rds")



####male
data.input.male <- GetAssayData(male, assay = "RNA", slot = "data")
labels <- Idents(male)
identity.male <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.male <- createCellChat(object = data.input.male)
cellchat.seu.male <- addMeta(cellchat.seu.male, meta = identity.male, 
                             meta.name = "labels")
cellchat.seu.male <- setIdent(cellchat.seu.male, ident.use = "labels") 

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.male@DB <- CellChatDB.use

cellchat.seu.male <- subsetData(cellchat.seu.male)
cellchat.seu.male <- identifyOverExpressedGenes(cellchat.seu.male)
cellchat.seu.male <- identifyOverExpressedInteractions(cellchat.seu.male)
cellchat.seu.male <- projectData(cellchat.seu.male, PPI.human)
cellchat.seu.male <- computeCommunProb(cellchat.seu.male)
cellchat.seu.male <- filterCommunication(cellchat.seu.male, min.cells = 10)
df.net.male <- subsetCommunication(cellchat.seu.male, slot.name = "net")
write.csv(df.net.male, "net.LR.male.csv")

cellchat.seu.male <- computeCommunProbPathway(cellchat.seu.male)
df.netp.male <- subsetCommunication(cellchat.seu.male, slot.name = "netP")
write.csv(df.netp.male, "net.patway.male.csv")

cellchat.seu.male <- aggregateNet(cellchat.seu.male)
cellchat.seu.male <- netAnalysis_computeCentrality(cellchat.seu.male)
saveRDS(cellchat.seu.male, file = "cellchat.male.rds")




object.list <- list(female = cellchat.female, 
                    male = cellchat.male)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))



pdf(file = "female.pdf", height = 10, width = 10)
groupSize_female <- as.numeric(table(cellchat.female@idents))
netVisual_circle(cellchat.female@net$weight, vertex.weight = groupSize_female, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength - female", 
                 arrow.width = 0.7, arrow.size = 0.7)
dev.off()



pdf(file = "male.pdf", height = 10, width = 10)
groupSize_male <- as.numeric(table(cellchat.male@idents))
netVisual_circle(cellchat.male@net$weight, vertex.weight = groupSize_male, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength - male", 
                 arrow.width = 0.7, arrow.size = 0.7)
dev.off()


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + 
    colSums(x@net$count)-diag(x@net$count)})

# control the dot size in the different datasets
weight.MinMax <- c(min(num.link), max(num.link)) 
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], 
                                               title = names(object.list)[i], 
                                               weight.MinMax = weight.MinMax)
}

pdf(file = "42整体互作.pdf", width = 10, height = 5)
patchwork::wrap_plots(plots = gg, nrow = 1)
dev.off()



cols_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
groupSize <- as.numeric(table(cellchat.male@idents))
mat <- cellchat.male@net$count
pdf(file = "51male.pdf", height = 15, width = 15)
par(mfrow = c(3,3), xpd=TRUE, mar=c(2,2,2,2))
for (i in 1:nrow(mat)){
  mat2 <- matrix (0, nrow = nrow (mat), ncol = ncol (mat), dimnames = dimnames (mat))
  mat2 [i, ] <- mat[i, ] 
  mat2 <- as.data.frame (mat2) 
  data <- data.frame(mat2 [i, ]) 
  data <- select(data, i, everything())
  datal <- data.frame (mat2[-i, ])
  datal <- select(datal, i, everything ())
  data2 <- rbind(data, datal)
  data2 <- as.matrix(data2)
  netVisual_circle(data2,
                   vertex.weight = groupSize,
                   weight.scale = T,
                   edge.weight.max = max(mat),
                   title.name = rownames(data2)[1],
                   color.use = cols_clusters,
                   label.edge = T,
                   edge.label.cex = 1.2,
                   vertex.label.cex = 1.2)
}
dev.off()




cols_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
groupSize <- as.numeric(table(cellchat.female@idents))
mat <- cellchat.female@net$count
pdf(file = "51female.pdf", height = 15, width = 15)
par(mfrow = c(3,3), xpd=TRUE, mar=c(2,2,2,2))
for (i in 1:nrow(mat)){
  mat2 <- matrix (0, nrow = nrow (mat), ncol = ncol (mat), dimnames = dimnames (mat))
  mat2 [i, ] <- mat[i, ] 
  mat2 <- as.data.frame (mat2) 
  data <- data.frame(mat2 [i, ]) 
  data <- select(data, i, everything())
  datal <- data.frame (mat2[-i, ])
  datal <- select(datal, i, everything ())
  data2 <- rbind(data, datal)
  data2 <- as.matrix(data2)
  netVisual_circle(data2,
                   vertex.weight = groupSize,
                   weight.scale = T,
                   edge.weight.max = max(mat),
                   title.name = rownames(data2)[1],
                   color.use = cols_clusters,
                   label.edge = T,
                   edge.label.cex = 1.2,
                   vertex.label.cex = 1.2)
}
dev.off()


pdf(file = "内皮调控成其他细胞.pdf", height = 10, width = 7)
netVisual_bubble(cellchat, sources.use = c(5), targets.use = c(1, 2, 3, 4, 5, 6, 
                                                               7, 8, 9, 10, 11),  
                 comparison = c(1, 2), angle.x = 45)
dev.off()


pathways.show <- "CXCL"
pdf(file = "female-CXCL.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat.female, signaling = pathways.show, layout = "chord")
dev.off()

pdf(file = "male-CXCL.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat.male, signaling = pathways.show, layout = "chord")
dev.off()



pathways.show <- "VEGF"
pdf(file = "female-VEGF.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat.female, signaling = pathways.show, layout = "chord")
dev.off()

pdf(file = "male-VEGF.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat.male, signaling = pathways.show, layout = "chord")
dev.off()



library(NMF)
library(ggalluvial)
pdf(file = "female-out.pdf", width = 10, height = 5)
selectK(cellchat.female, pattern = "outgoing")
dev.off()


pdf(file = "female-in.pdf", width = 10, height = 5)
selectK(cellchat.female, pattern = "incoming")
dev.off()


pdf(file = "male-in.pdf", width = 10, height = 5)
selectK(cellchat.male, pattern = "incoming")
dev.off()

pdf(file = "male-out.pdf", width = 10, height = 5)
selectK(cellchat.male, pattern = "outgoing")
dev.off()


pdf(file = "female-pattern-out.pdf", width = 10, height = 8)
identifyCommunicationPatterns(object = cellchat.female, pattern = "outgoing", k = 2, 
                              height = 10)
dev.off()


# river plot
pdf(file = "female.out.river.pdf", width = 10, height = 9)
netAnalysis_river(cellchat.female, pattern = "outgoing")
dev.off()


# dot plot
pdf(file = "female.out.dot.pdf", width = 10, height = 9)
netAnalysis_dot(cellchat.female, pattern = "outgoing")
dev.off()

pdf(file = "female-pattern-in.pdf", width = 10, height = 8)
identifyCommunicationPatterns(object = cellchat.female, 
                              pattern = "incoming", k = 2, height = 10)
dev.off()


cellchat.female <- identifyCommunicationPatterns(object = cellchat.female, pattern = "incoming", k = 2, height = 10)
pdf(file = "female.in.river.pdf", width = 10, height = 9)
netAnalysis_river(cellchat.female, pattern = "incoming")
dev.off()

pdf(file = "female.in.dot.pdf", width = 10, height = 9)
netAnalysis_dot(cellchat.female, pattern = "incoming")
dev.off()


cellchat.male <- identifyCommunicationPatterns(object = cellchat.male, pattern = "outgoing", k = 3)
pdf(file = "male-pattern-out.pdf", width = 10, height = 8)
identifyCommunicationPatterns(object = cellchat.male, pattern = "outgoing", k = 3, height = 10)
dev.off()


pdf(file = "male.out.river.pdf", width = 10, height = 9)
netAnalysis_river(cellchat.male, pattern = "outgoing")
dev.off()


pdf(file = "male.out.dot.pdf", width = 10, height = 9)
netAnalysis_dot(cellchat.male, pattern = "outgoing")
dev.off()


cellchat.male <- identifyCommunicationPatterns(object = cellchat.male, pattern = "incoming", k = 2)
pdf(file = "male-pattern-in.pdf", width = 10, height = 8)
identifyCommunicationPatterns(object = cellchat.male, pattern = "incoming", k = 2, height = 10)
dev.off()

pdf(file = "male.in.river.pdf", width = 10, height = 9)
netAnalysis_river(cellchat.male, pattern = "incoming")
dev.off()








































































































