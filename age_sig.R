library(Seurat)
library(CellChat)
table(seu$Age2)
library(CellChat)
library(Seurat)
table(seu$celltype)

options(stringsAsFactors = FALSE)
options(futrue.globlas.Maxsize=2000*1024**3)
conflicts_prefer(SeuratObject::saveRDS)



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

DefaultAssay(seu) <- "RNA"
seu$seurat_clusters <- NULL
seu$SCT_snn_res.0.45 <- NULL
seu$SCT_snn_res.0.4 <- NULL
seu$SCT_snn_res.0.35 <- NULL
seu$SCT_snn_res.0.3 <- NULL
seu$SCT_snn_res.0.25 <- NULL
seu$SCT_snn_res.0.2 <- NULL
seu$SCT_snn_res.0.15 <- NULL
seu$SCT_snn_res.0.1 <- NULL
seu$SCT_snn_res.0.05 <- NULL
seu$SCT_snn_res.0.05 <- NULL
seu$SCT_snn_res.0 <- NULL
seu$nCount_SCT <- NULL
seu$nFeature_SCT <- NULL
seu$percent.mt <- NULL
seu$G2M.Score <- NULL
seu$S.Score <- NULL
seu$Phase <- NULL
seu@assays$SCT <- NULL
seu@reductions$pca <- NULL
seu@reductions$harmony <- NULL
seu@reductions$umap <- NULL

table(Y$celltype1)
table(M$celltype1)
table(O$celltype1)



Y <- subset(seu, Age2 == "Y")
M <- subset(seu, Age2 == "M")
O <- subset(seu, Age2 == "O")


####Y
###Y
data.input.Y <- GetAssayData(Y, assay = "RNA", slot = "data")
labels <- Idents(Y)
identity.Y <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.Y <- createCellChat(object = data.input.Y)
cellchat.seu.Y <- addMeta(cellchat.seu.Y, meta = identity.Y, 
                          meta.name = "labels")
cellchat.seu.Y <- setIdent(cellchat.seu.Y, ident.use = "labels") 

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.Y@DB <- CellChatDB.use

cellchat.seu.Y <- subsetData(cellchat.seu.Y)
cellchat.seu.Y <- identifyOverExpressedGenes(cellchat.seu.Y)
cellchat.seu.Y <- identifyOverExpressedInteractions(cellchat.seu.Y)
cellchat.seu.Y <- projectData(cellchat.seu.Y, PPI.human)
cellchat.seu.Y <- computeCommunProb(cellchat.seu.Y)
cellchat.seu.Y <- filterCommunication(cellchat.seu.Y, min.cells = 10)
df.net.Y <- subsetCommunication(cellchat.seu.Y)
write.csv(df.net.Y, "net.LR.Y.csv")

cellchat.seu.Y <- computeCommunProbPathway(cellchat.seu.Y)
df.netp.Y <- subsetCommunication(cellchat.seu.Y, slot.name = "netP")
write.csv(df.netp.Y, "net.patway.Y.csv")

cellchat.seu.Y <- aggregateNet(cellchat.seu.Y)
cellchat.seu.Y <- netAnalysis_computeCentrality(cellchat.seu.Y)
saveRDS(cellchat.seu.Y, file = "cellchat.Y.rds")


####M
###M
data.input.M <- GetAssayData(M, assay = "RNA", slot = "data")
labels <- Idents(M)
identity.M <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.M <- createCellChat(object = data.input.M)
cellchat.seu.M <- addMeta(cellchat.seu.M, meta = identity.M, 
                          meta.name = "labels")
cellchat.seu.M <- setIdent(cellchat.seu.M, ident.use = "labels") 

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.M@DB <- CellChatDB.use

cellchat.seu.M <- subsetData(cellchat.seu.M)
cellchat.seu.M <- identifyOverExpressedGenes(cellchat.seu.M)
cellchat.seu.M <- identifyOverExpressedInteractions(cellchat.seu.M)
cellchat.seu.M <- projectData(cellchat.seu.M, PPI.human)
cellchat.seu.M <- computeCommunProb(cellchat.seu.M)
cellchat.seu.M <- filterCommunication(cellchat.seu.M, min.cells = 10)
df.net.M <- subsetCommunication(cellchat.seu.M)
write.csv(df.net.M, "net.LR.M.csv")


cellchat.seu.M <- computeCommunProbPathway(cellchat.seu.M)
df.netp.M <- subsetCommunication(cellchat.seu.M, slot.name = "netP")
write.csv(df.netp.M, "net.patway.M.csv")

cellchat.seu.M <- aggregateNet(cellchat.seu.M)
cellchat.seu.M <- netAnalysis_computeCentrality(cellchat.seu.M)
saveRDS(cellchat.seu.M, file = "cellchat.M.rds")


#####O
###O
data.input.O <- GetAssayData(O, assay = "RNA", slot = "data")
labels <- Idents(O)
identity.O <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.O <- createCellChat(object = data.input.O)
cellchat.seu.O <- addMeta(cellchat.seu.O, meta = identity.O, 
                          meta.name = "labels")
cellchat.seu.O <- setIdent(cellchat.seu.O, ident.use = "labels") 

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.O@DB <- CellChatDB.use

cellchat.seu.O <- subsetData(cellchat.seu.O)
cellchat.seu.O <- identifyOverExpressedGenes(cellchat.seu.O)
cellchat.seu.O <- identifyOverExpressedInteractions(cellchat.seu.O)
cellchat.seu.O <- projectData(cellchat.seu.O, PPI.human)
cellchat.seu.O <- computeCommunProb(cellchat.seu.O)
cellchat.seu.O <- filterCommunication(cellchat.seu.O, min.cells = 10)
df.net.O <- subsetCommunication(cellchat.seu.O)
write.csv(df.net.O, "net.LR.O.csv")

cellchat.seu.O <- computeCommunProbPathway(cellchat.seu.O)
df.netp.O <- subsetCommunication(cellchat.seu.O, slot.name = "netP")
write.csv(df.netp.O, "net.patway.O.csv")

cellchat.seu.O <- aggregateNet(cellchat.seu.O)
cellchat.seu.O <- netAnalysis_computeCentrality(cellchat.seu.O)
saveRDS(cellchat.seu.O, file = "cellchat.O.rds")



object.list <- list(Y = cellchat.seu.Y, 
                    M = cellchat.seu.M,
                    O = cellchat.seu.O)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


pdf(file = "Y.pdf", height = 10, width = 10)
groupSize_Y <- as.numeric(table(cellchat.seu.Y@idents))
netVisual_circle(cellchat.seu.Y@net$weight, vertex.weight = groupSize_Y, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength - Y", 
                 arrow.width = 0.7, arrow.size = 0.7)
dev.off()



pdf(file = "M.pdf", height = 10, width = 10)
groupSize_M <- as.numeric(table(cellchat.seu.M@idents))
netVisual_circle(cellchat.seu.M@net$weight, vertex.weight = groupSize_M, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength - M", 
                 arrow.width = 0.7, arrow.size = 0.7)
dev.off()


pdf(file = "O.pdf", height = 10, width = 10)
groupSize_O <- as.numeric(table(cellchat.seu.O@idents))
netVisual_circle(cellchat.seu.O@net$weight, vertex.weight = groupSize_O, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength - O", 
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

pdf(file = "整体互作.pdf", width = 15, height = 5)
patchwork::wrap_plots(plots = gg, nrow = 1)
dev.off()


cols_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
groupSize <- as.numeric(table(cellchat.seu.Y@idents))
mat <- cellchat.seu.Y@net$count
pdf(file = "Y-vec.pdf", height = 15, width = 15)
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
groupSize <- as.numeric(table(cellchat.seu.M@idents))
mat <- cellchat.seu.M@net$count
pdf(file = "M-vec.pdf", height = 15, width = 15)
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
groupSize <- as.numeric(table(cellchat.seu.O@idents))
mat <- cellchat.seu.O@net$count
pdf(file = "O-vec.pdf", height = 15, width = 15)
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


pdf(file = "内皮调控成其他细胞.pdf", height = 10, width = 10)
netVisual_bubble(cellchat, sources.use = c(5), targets.use = c(1, 2, 3, 4, 5, 6, 
                                                               7, 8, 9, 10, 11),  
                 comparison = c(1, 2, 3), angle.x = 45)
dev.off()






































