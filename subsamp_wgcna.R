
# 设置
args <- commandArgs(TRUE)
norm_file_dir <- "/cbcb/lab/smount/ZCL/bioconductor_scripts/counts/subsets/"

# 获取命令行参数
if (length(args) < 7) {
  message("参数不足")
  quit(save = "no")
} else {
  softPower <- as.numeric(args[1])  # 平滑参数
  minModuleSize <- as.numeric(args[2])  # 最小模块大小
  merge_eigengenes <- as.numeric(args[3])  # 是否合并特征基因（1 = 是，其他 = 否）
  gene_sub_list_name <- as.character(args[4])  # 基因子集文件名
  norm_file <- as.character(args[5])  # 规范化文件名
  save_dir <- as.character(args[6])  # 保存目录
  runNum <- as.character(args[7])  # 运行编号，用于跟踪进程成功情况
}

# 规范化文件路径
filename <- paste(norm_file_dir, norm_file, sep = "")

# 创建失败标记文件
file.create(paste("/cbcb/project-scratch/ZCL/wgcna/consensus/", save_dir, "/failed/", runNum, sep = ""))

# 设置工作目录
setwd(paste("/cbcb/project-scratch/ZCL/wgcna/consensus/", save_dir, sep = ""))

# 重定向输出
# sink(file = paste(getwd(), "/output/", runNum, ".txt", sep = ""))

# 打印平滑参数
# print(paste("softPower: ", , sep = ""))

# 加载库
suppressMessages(library("WGCNA"))
suppressMessages(library("methods"))
suppressMessages(library("edgeR"))
# suppressMessages(library(gplots))
suppressMessages(library(RColorBrewer))
nthr <- enableWGCNAThreads()  # 启用多线程
print(nthr)

print("导入数据")
selection <- read.csv(filename, header = TRUE)

# 设置行名
rownames(selection) <- selection[,"X"]
# 删除无用列
selection[,"X"] <- NULL
selection[,"X.1"] <- NULL

print(head(selection))

####
# 打开基因子集文件并子集名称
####

sub_genes <- read.delim(paste("sub_genes/", gene_sub_list_name, sep = ""), sep = "\n", header = FALSE)
print(sub_genes$V1)
selection <- selection[as.character(sub_genes$V1), ]
print(head(selection))

print(dim(selection))
# 转置矩阵并保留数值
transpose <- t(selection)
# selection <- NULL

# "gene35204-v1.0-hybrid"

# 移除非方差列
print("移除方差很小的列")
transpose <- transpose[, apply(transpose, 2, var, na.rm = TRUE) > 0.05]
print(dim(transpose))

print("生成邻接矩阵")
adjacency_res <- adjacency(transpose, power = softPower, type = "signed")

dim(adjacency_res)

distanceMat <- 1 - adjacency_res

print("构建基因树")
# 调用层次聚类函数
geneTree <- hclust(as.dist(distanceMat), method = "average")

# 使用动态树切割模块
print("切割模块")
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = distanceMat, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize, verbose = 4)

table(dynamicMods)

# 将数字标签转换为颜色
dynamicColors <- labels2colors(dynamicMods)
l <- length(unique(dynamicColors))
l
table(dynamicColors)

if (l > 200) {
  file.remove(paste("/cbcb/project-scratch/ZCL/wgcna/consensus/", save_dir, "/failed/", runNum, sep = ""))
  print("模块过多")
  file.create(paste("/cbcb/project-scratch/ZCL/wgcna/consensus/", save_dir, "/2many/", runNum, sep = ""))
  quit(save = "no")
}

# 写出簇
gene_names <- rownames(adjacency_res)
clusters <- as.data.frame(cbind(genes = gene_names, group = dynamicColors))
clusters <- clusters[order(clusters$group), ]

if (merge_eigengenes == 0) {
  write.csv(clusters, file = paste("/cbcb/project2-scratch/ZCL/", save_dir, "/clusters/", runNum, ".csv", sep = ""), row.names = FALSE)
}

moduleColors <- dynamicColors
# 如果需要基于特征基因合并簇
if (merge_eigengenes == 1) {
  # 计算特征基因
  MEList <- moduleEigengenes(transpose, colors = dynamicColors)
  MEs <- MEList$eigengenes
  # 计算模块特征基因的不相似性
  MEDiss <- 1 - cor(MEs)
  # 聚类模块特征基因
  METree <- hclust(as.dist(MEDiss), method = "average")

  MEDissThres <- 0.25

  # 调用自动合并函数
  merge <- mergeCloseModules(transpose, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # 合并后的模块颜色
  mergedColors <- merge$colors
  moduleColors <- mergedColors
  print(table(moduleColors))
  print(length(unique(moduleColors)))
  # 合并模块的特征基因
  mergedMEs <- merge$newMEs

  # 写出合并后的簇
  gene_names <- rownames(adjacency_res)
  mclusters <- as.data.frame(cbind(genes = gene_names, group = mergedColors))
  mclusters <- mclusters[order(mclusters$group), ]
  write.csv(mclusters, file = paste("/cbcb/project2-scratch/ZCL/", save_dir, "/clusters/", runNum, ".csv", sep = ""), row.names = FALSE)
}

print("完成")

# 停止输出重定向
# sink()

file.remove(paste("/cbcb/project-scratch/ZCL/wgcna/consensus/", save_dir, "/failed/", runNum, sep = ""))
file.create(paste("/cbcb/project-scratch/ZCL/wgcna/consensus/", save_dir, "/success/", runNum, sep = ""))
