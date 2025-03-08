"""
# 输入参数说明
# - 实验目录名称
# - 基因名称列表
# - 输入计数矩阵名称
# - 采样次数
# - 布尔值 ("rand" = T, 其他值 = F) 描述是否使用随机参数
#   否则使用默认的硬编码参数
output: experiment directory is created in /cbcb/project-scratch/ZCL/wgcna/
		subdirectories sub_genes/, clusters/, ouput/, failed/, success/, 2many/, and bash/ are created
		file containing the number of times a gene was chosen and the number of subsamplings that were taken

	programs takes a subsample of gene names for the network and some parameters for network creation and produces bash file to run the R script subsamp_wgcna.R in the bash/ subdirectory of the experiment directory
"""

import sys
import os
import random
import subprocess

# 如果命令行参数中包含 "-ls"，则列出所有目录
if "-ls" in sys.argv:
    subprocess.call(["ls", "-d", "*/"])

# 获取命令行参数
exp_dir = sys.argv[1]  # 实验目录名称
count_mat = sys.argv[2]  # 输入的计数矩阵名称
n = int(sys.argv[3])  # 采样次数
rand = sys.argv[4]  # 是否使用随机参数

# 设置采样率
sub_samp_rate = 0.8

print("making directories")

# 尝试创建实验目录及其子目录
try:
    os.mkdir("/cbcb/project2-scratch/ZCL/" + exp_dir)  # 创建实验目录
    os.mkdir("/cbcb/project2-scratch/ZCL/" + exp_dir + "/clusters")  # 创建 clusters 子目录
    os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir)  # 创建 consensus 子目录
    os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/ouput")  # 创建 ouput 子目录
    os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/failed")  # 创建 failed 子目录
    os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/success")  # 创建 success 子目录
    os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/2many")  # 创建 2many 子目录
    os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/bash")  # 创建 bash 子目录
    os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/sub_genes")  # 创建 sub_genes 子目录
    os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/clusters")  # 创建 clusters 子目录
    os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/config")  # 创建 config 子目录
    os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/adjmat")  # 创建 adjmat 子目录
    os.mkdir("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/indmat")  # 创建 indmat 子目录

except:
    print("experiment name already used")  # 如果目录已存在，打印错误信息
    exit(1)  # 退出程序

# 初始化基因使用情况字典和基因列表
gene_used = {}
gene_list = []

# 打开基因列表文件
genes = open("/cbcb/lab/smount/ZCL/gene_list.txt")

print("reading gene file")
# 读取基因列表文件中的每一行
for g in genes:
    gene_list.append(g.strip("\n"))  # 去除换行符并添加到基因列表
    gene_used.update({g.strip("\n"): 0})  # 初始化基因使用次数为0

gene_num = len(gene_list)  # 基因总数
select_frac = int(gene_num * 0.8)  # 选择的基因数量

# 定义网络构建参数
powers = [1, 2, 4, 8, 12, 16]
minModuleSize = [40, 60, 90, 120, 150, 180, 210]
if count_mat[0] == "l":  # 如果计数矩阵名称以 "l" 开头，调整 minModuleSize
    minModuleSize = [90, 120, 150, 180, 210]
merge_eigengene = [0, 1]

# 打开配置文件
cfg = open("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/cluster_config.txt", 'w')
cfg.write("power\tminModSize\tmerge\n")  # 写入配置文件头部

print("generating gene combinations")
# 生成基因组合
for x in range(0, n):
    samp = random.sample(gene_list, select_frac)  # 随机选择基因
    print("sample size: " + str(len(samp)))  # 打印样本大小

    # 创建样本文件
    f = open("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/sub_genes/" + "sample" + str(x) + ".txt", 'w')
    for s in samp:
        f.write(s + "\n")  # 写入样本基因
        gene_used[s] += 1  # 更新基因使用次数
    f.close()

    # 根据是否使用随机参数选择网络构建参数
    if rand == "rand":
        random.shuffle(powers)  # 随机打乱 powers
        random.shuffle(minModuleSize)  # 随机打乱 minModuleSize
        random.shuffle(merge_eigengene)  # 随机打乱 merge_eigengene
        p = powers[0]  # 选择第一个 powers
        minM = minModuleSize[0]  # 选择第一个 minModuleSize
        m = merge_eigengene[0]  # 选择第一个 merge_eigengene
    else:
        p = 2  # 使用默认 power
        minM = 90  # 使用默认 minModuleSize
        m = 0  # 使用默认 merge

    # 写入配置文件
    cfg.write(str(p) + "\t" + str(minM) + "\t" + str(m) + "\n")

    # 创建运行脚本文件
    run = open("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/bash/run" + str(x) + ".sh", 'w')
    run.write("#PBS -q throughput\n#PBS -l mem=36GB,walltime=12:00:00,ncpus=2\n")  # 写入 PBS 资源请求
    run.write("/cbcb/lab/smount/programs/R-3.1.2/bin/Rscript /cbcb/lab/smount/ZCL/bioconductor_scripts/subsamp_wgcna.R " + str(p) + " " + str(minM) + " " + str(m) + " " + "sample" + str(x) + ".txt" + " " + count_mat + " " + exp_dir + " " + str(x) + "\n")  # 写入 R 脚本调用
    run.close()

cfg.close()

print("writing gene sample rates")
# 创建基因采样率文件
gene_file = open("/cbcb/project-scratch/ZCL/wgcna/consensus/" + exp_dir + "/gene_sample_rate.csv", 'w')
gene_file.write("gene\tselected\ttotal_sample\n")  # 写入文件头部
for k in gene_used.keys():
    gene_file.write(k + "\t" + str(gene_used[k]) + "\t" + str(n) + "\n")  # 写入每个基因的采样情况
gene_file.close()