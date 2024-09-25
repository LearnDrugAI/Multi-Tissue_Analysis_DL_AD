# 加载必要的包
library(data.table)
library(UpSetR)

# 读取数据
data_father <- fread("CMAP_combined_with_qvalue_zscore_GSCT90042665.csv")
data_mother <- fread("CMAP_combined_with_qvalue_zscore_GSCT90042678.csv")
data_siblings <- fread("CMAP_combined_with_qvalue_zscore_GSCT90042691.csv")

# 获取药物和组织的交集
shared_drugs <- Reduce(intersect, list(data_father$Drug, data_mother$Drug, data_siblings$Drug))
shared_tissues <- Reduce(intersect, list(data_father$tissue, data_mother$tissue, data_siblings$tissue))

# 确保 Siblings_Tissues 的数据即使数量少也显示
listInput <- list(
  Father_Drugs = data_father$Drug,
  Mother_Drugs = data_mother$Drug,
  Siblings_Drugs = data_siblings$Drug,
  Father_Tissues = data_father$tissue,
  Mother_Tissues = data_mother$tissue,
  Siblings_Tissues = data_siblings$tissue
)

# 确保关闭所有现有的图形设备以防止生成空白页
graphics.off()

# 打开PDF设备
pdf("UpSet_Drug_Tissue_Intersection_with_Fixed_Counts.pdf", width = 10, height = 6)

# 创建UpSet图，确保显示药物和组织交集的准确数量，并显示所有集合
upset(
  fromList(listInput),
  order.by = "freq",
  main.bar.color = "steelblue",
  sets.bar.color = "darkgreen",
  text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5),
  keep.order = TRUE,
  sets.x.label = "Set Size",   # 设置x轴标签为Set Size
  show.numbers = "yes",        # 显示交集的数量
  set_size.show = TRUE,        # 确保显示集合大小
  set_size.scale_max = 2000,   # 设置集合大小刻度的最大值
  mb.ratio = c(0.6, 0.4),      # 设置主条形图和集合条形图的比例
  nsets = length(listInput)    # 强制显示所有集合，即使数据量少
)

# 关闭PDF设备
dev.off()

# 确保关闭所有图形设备
graphics.off()
