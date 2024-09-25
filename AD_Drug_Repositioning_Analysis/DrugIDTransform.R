#***********************************
# Convert the ID of the drug
#***********************************
# 加载必要的库
library(dplyr)
library(data.table)
library(stringr)

# 设置工作目录到包含csv文件的文件夹
setwd("~/mygenedata/GCST90042665")

# 获取所有文件名
file_names <- list.files(pattern = "^Results_5MethodAvg_GCST90042665_.*\\.csv$")

# 提取"Results_5MethodAvg_GCST90042678_"到".csv"之间的部分
extracted_names <- sub("^Results_5MethodAvg_GCST90042665_(.*)\\.csv$", "\\1", file_names)

# 设置工作目录到主目录
setwd("~/mygenedata")

# 设置输入和输出目录
input_dir <- "GCST90042665"
output_dir <- "GCST90042665"

# 创建输出目录如果不存在
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 定义组织列表
tissues <- extracted_names

# 获取组织的数量
no_tissue <- length(tissues)

# 遍历每个组织
for (i in 1:no_tissue) {
  # 构建文件名
  input_file <- file.path(input_dir, paste0("Results_5MethodAvg_GCST90042665_", tissues[i], ".csv"))
  output_file <- file.path(output_dir, paste0("results_5MethodAvg_GCST90042665_", tissues[i], ".csv"))
  
  # 尝试读取输入文件，如果出错则继续下一个
  try({
    # 读取结果文件
    reDrug <- fread(input_file, header = TRUE)
    reDrug <- reDrug[, -1]  # 删除第一列
    
    # 读取SigInfo文件
    SigInfo <- fread('GSE70138_Broad_LINCS_sig_info.txt')
    
    # 修改 reDrug 中 Drug 列的字符串
    reDrug <- reDrug %>%
      mutate(Drug = str_replace(Drug, "\\.(?=.{3}$)", ":"))
    
    # 创建SigInfo的子集
    SigInfo_subset <- SigInfo %>%
      dplyr::select(sig_id, pert_iname)
    
    # 匹配sig_id并替换Drug列为pert_iname
    result <- reDrug %>%
      left_join(SigInfo_subset, by = c("Drug" = "sig_id")) %>%
      dplyr::select(pert_iname, everything(), -Drug) %>%
      rename(Drug = pert_iname)
    
    # 写入结果到输出文件
    fwrite(result, output_file, row.names = FALSE)
    
    # 打印结果
    print(result)
  }, silent = TRUE)
}
