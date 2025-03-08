

#####
##### 这个脚本用于创建由 subsamp_wgcna.R 生成的簇的 1,0 指示矩阵
#####

library("data.table")
library(pryr)

# 获取命令行参数
args <- commandArgs(TRUE)

exp_name <- as.character(args[1])  # 实验名称
start <- as.numeric(args[2])      # 起始编号
stop <- as.numeric(args[3])       # 结束编号

## 从基因列表创建矩阵
gene_f <- read.delim("/cbcb/lab/smount/ZCL/gene_list.txt", sep = "\t", header = FALSE)

genes <- as.character(gene_f$V1)  # 基因名称列表

#############################################################
### 构造主矩阵作为 data.table
#############################################################

# 创建矩阵
print("make matrix")
master_mat <- matrix(0, nrow = length(genes), ncol = length(genes))

# 转换为数据框
print("make data frame")
master_df <- as.data.frame(master_mat)
rm(master_mat)

# 设置列名
colnames(master_df) <- genes
master_df$g <- genes

print("make data table")
# 转换为数据表
master_dt <- data.table(master_df)
rm(master_df)

setkey(master_dt, g)  # 设置键

master_dt[, g := NULL]  # 删除列 g

#############################################################
### 构造矩阵作为 data.table
#############################################################

# 创建矩阵
print("make matrix")
con_mat <- matrix(0, nrow = length(genes), ncol = length(genes))
print(dim(con_mat))

# 转换为数据框
print("make data frame")
con_df <- as.data.frame(con_mat)
rm(con_mat)

# 设置列名
colnames(con_df) <- genes

print("make data table")
# 转换为数据表
con_dt <- data.table(con_df)
rm(con_df)

# setkey(con_dt, g)

#############################################################
#############################################################
#############################################################

# 循环处理从 start 到 stop 的每个样本
for(i in seq(start, stop)){
    print(i)

    con_dt[,] = 0  # 重置矩阵
    con_dt[, g := genes]  # 添加基因列
    setkey(con_dt, g)  # 设置键

    exp_dir <- "/cbcb/project-scratch/ZCL/wgcna/consensus/"  # 实验目录

    # 打开簇列表文件
    group_name <- paste(exp_dir, exp_name, "/sub_genes/sample", i, ".txt", sep = "")
    group <- as.character(read.csv(group_name, header = FALSE)$V1)  # 读取基因列表

    print("filling matrix")
    con_dt[group, (group) := 1]  # 填充矩阵
    print("done filling matrix")

    print(mem_used())  # 打印内存使用情况

    con_dt[, g := NULL]  # 删除列 g
    master_dt <- master_dt + con_dt  # 累加到主矩阵
}

# 写入主矩阵到文件
fwrite(master_dt, file = paste(exp_dir, exp_name, "/indmat/", start, "_", stop, ".csv", sep = ""))


