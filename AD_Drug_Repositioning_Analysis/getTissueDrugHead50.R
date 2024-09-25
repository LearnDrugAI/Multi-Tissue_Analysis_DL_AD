# 加载必要的库
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(WGCNA)
library(qvalue)
library(grid)
library(ggtext)
# 设置输入文件目录
input_dir <- "GSCT90042691"

# 获取文件名列表
file_names <- list.files(path = "GSCT90042691/transIDtoName", pattern = "^results_5MethodAvg_GCST90042691_.*\\.csv$", full.names = TRUE)

# 提取组织名（tissues）
extracted_names <- sub("^GSCT90042691/transIDtoName/results_5MethodAvg_GCST90042691_(.*)\\.csv$", "\\1", file_names)
tissues <- extracted_names
no_tissue <- length(tissues)

# 初始化一个空列表用于存储所有组织的数据
CMAP_combined <- list()

# 读取并合并所有组织的数据，保留 Drug、AvgRank 和 perm.p 列
for (i in 1:no_tissue) {
  CMAP <- read.csv(file_names[i])  # 直接使用file_names中的路径读取文件
  if (nrow(CMAP) > 0) {  # 确保文件不为空
    CMAP <- CMAP[, c("Drug", "AvgRank", "perm.p")]  # 保留 Drug、AvgRank 和 perm.p 列
    CMAP$qvalue <- qvalue(CMAP$perm.p)$qvalues  # 计算q-value
    CMAP$zscore <- (CMAP$perm.p - mean(CMAP$perm.p, na.rm = TRUE)) / sd(CMAP$perm.p, na.rm = TRUE)  # 计算z-score
    CMAP$tissue <- tissues[i]  # 添加一个列标识当前的组织
    CMAP_combined[[i]] <- CMAP  # 将数据存入列表
  }
}
# 将列表转换为数据框
CMAP_combined <- do.call(rbind, CMAP_combined)
# 设置保存文件的路径
output_file <- "CMAP_combined_with_qvalue_zscore_GSCT90042691.csv"

# 将数据保存到文件中
write.csv(CMAP_combined, file = output_file, row.names = FALSE)

# 检查是否成功读取了数据
if (length(CMAP_combined) == 0) {
  stop("未能读取任何数据，请检查文件路径和文件内容。")
}


# # 选择前50个qvalue最小的药物
# CMAP_qvalue_sorted <- CMAP_combined[order(CMAP_combined$qvalue), ]

# 按照 qvalue 和 perm.p 进行排序
CMAP_qvalue_sorted <- CMAP_combined[order(CMAP_combined$qvalue, CMAP_combined$perm.p), ]
top50_drugs_qvalue <- CMAP_qvalue_sorted[1:50]
# # 初始化一个空向量来存储唯一的药物
# unique_drugs <- c()
# 
# # 初始化计数器来记录处理到第几行时有50个独特的药物
# cutoff_index <- 0
# 
# # 遍历排序后的数据框，直到我们找到50个不同的药物
# for (i in 1:nrow(CMAP_qvalue_sorted)) {
#   drug <- CMAP_qvalue_sorted$Drug[i]
#   if (!(drug %in% unique_drugs)) {
#     unique_drugs <- c(unique_drugs, drug)
#   }
#   if (length(unique_drugs) == 50) {
#     cutoff_index <- i
#     break
#   }
# }
# 
# # 现在我们知道到第 cutoff_index 行为止包含50个独特的药物
# # 所以直接选择 CMAP_qvalue_sorted 的前 cutoff_index 行
# top50_drugs_qvalue <- CMAP_qvalue_sorted[1:cutoff_index, ]
## 从原始数据中过滤出这些排名前50的药物及其对应的组织
#CMAP_filtered_qvalue <- subset(CMAP_qvalue_sorted, Drug %in% top50_drugs_qvalue$Drug)

# 转换为长格式（tidy data），只保留 Drug 和 qvalue 列
CMAP_melted_qvalue <- melt(top50_drugs_qvalue, id.vars = c("Drug", "tissue"), measure.vars = "qvalue")


# 调色板设置
#my_palette <- colorRampPalette(c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B"), space = "Lab")
#myPalette2way <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")), space = "Lab")
my_palette <- colorRampPalette(rev(c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B")), space = "Lab")
#my_palette <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)

# 自定义PDF函数：绘制基于某个统计量
# save_heatmap_to_pdf <- function(CMAP_melted, fill_var, pdf_file, base_width = 10, base_height = 8) {
#   n_drugs <- length(unique(CMAP_melted$Drug))
#   n_tissues <- length(unique(CMAP_melted$tissue))
#  
#   width <- base_width + n_drugs * 0.5
#   height <- base_height + n_tissues * 0.5
#   
#   pdf(pdf_file, width = width, height = height)
#   
#   p <- ggplot(CMAP_melted, aes(x = Drug, y = tissue, fill = value)) +
#     geom_tile(color = "white") +
#     scale_fill_gradientn(colours = my_palette(100)) +
#     theme_minimal() +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1, size = 20),  # 增加X轴标签字体大小
#       axis.text.y = element_text(size = 20),  # 增加Y轴标签字体大小
#       axis.title.x = element_text(size = 22, face = "bold"),  # 增加X轴标题字体大小并加粗
#       axis.title.y = element_text(size = 22, face = "bold"),  # 增加Y轴标题字体大小并加粗
#       plot.title = element_text(size = 24, face = "bold", hjust = 0.5),  # 增加图标题字体大小并加粗
#       legend.title = element_text(size = 24),  # 增加图例标题字体大小
#       legend.text = element_text(size = 22)  # 增加图例标签字体大小
#     ) +
#     guides(fill = guide_colorbar(barheight = 15)) +  # 增加图例的高度
#     labs(x = "Drug", y = "Tissue", fill = fill_var, title = paste("Top 50 Repurposed Drugs for Illnesses of siblings:Alzheimer's disease/dementia Based on Q-value in Significant Tissues"))
#   
#   print(p)
#   dev.off()
# } 



# # 调用函数，生成基于 q-value 的热图并保存为 PDF
# save_heatmap_to_pdf(CMAP_melted_qvalue, "qvalue", pdf_file = "Tissue_Drug_Associations_Heatmap_top50_qvalue_2691.pdf")


# 生成基于 z-score 的热图并保存为 PDF
#save_heatmap_to_pdf(CMAP_melted_zscore, "zscore", pdf_file = "Tissue_Drug_Associations_Heatmap_top50_zscore_1866.pdf")
# 假设 top50_drugs_qvalue 是已经包含前50个药物的筛选结果

# 转换数据框，准备绘制分组条形图
# 计算每种药物的重复次数

top50_drugs_qvalue$Drug_count <- ave(top50_drugs_qvalue$Drug, top50_drugs_qvalue$Drug, FUN = length)
# 设置重复药物的颜色
highlight_color <- "#67001F"  
# 将 Drug_count 转换为因子，以便正确显示在图例中
top50_drugs_qvalue$Drug_count <- factor(top50_drugs_qvalue$Drug_count)
# 首先将数据从宽格式转换为长格式
# Convert p-values and q-values to -log10 scale
top50_drugs_qvalue$log10_perm.p <- -log10(top50_drugs_qvalue$perm.p)
top50_drugs_qvalue$log10_qvalue <- -log10(top50_drugs_qvalue$qvalue)
# 转换为长格式
CMAP_long_log <- melt(top50_drugs_qvalue, id.vars = c("Drug", "tissue", "Drug_count"), 
                      measure.vars = c("log10_perm.p", "log10_qvalue"))
# 找到 P-value 中的最高值和最低值
min_value <- min(top50_drugs_qvalue$log10_perm.p)
max_value <- max(top50_drugs_qvalue$log10_perm.p)

# 获取对应最小和最大 P-value 的药物名称
min_drug <- top50_drugs_qvalue$Drug[top50_drugs_qvalue$log10_perm.p == min_value]
max_drug <- top50_drugs_qvalue$Drug[top50_drugs_qvalue$log10_perm.p == max_value]

# 获取 Q-value 的 log10 转换值
qvalue_label <- round(top50_drugs_qvalue$log10_qvalue[1], 2)

# 创建条形图
p <- ggplot(CMAP_long_log, aes(x = reorder(Drug, value), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  scale_fill_manual(values = c("#92C5DE", "#41AB5D"), 
                    labels = c(expression(-log[10]("P-value")), 
                               expression(-log[10]("Q-value")))) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12, face = "bold", margin = margin(r = 20)),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  ) +
  labs(
    title = "Top 50 Repurposed Drugs for Illnesses of siblings:Alzheimer's disease/dementia  \nBased on Q−value and P-value in Stomach",
    x = "Drug",
    y = expression(-log[10]("Value"))
  ) +
  # 只标注 P-value 的最高和最低值
  geom_text(data = subset(CMAP_long_log, variable == "log10_perm.p" & 
                            (Drug == min_drug | Drug == max_drug)),
            aes(label = round(value, 2)), 
            position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
  # 添加一条 -log10(Q-value) 的虚线
  geom_hline(yintercept = qvalue_label, linetype = "dashed", color = "black", size = 1, alpha = 0.7) +
  
  # 自定义 y 轴刻度标签
  scale_y_continuous(
    breaks = c(seq(0, max_value, by = 1), qvalue_label),
    labels = function(x) {
      sapply(x, function(val) {
        if (val == qvalue_label) {
          # return(glue::glue("<b>{qvalue_label}</b>"))   # 使用HTML加粗qvalue_label
          return(NULL)
        } else {
          return(as.character(val))
        }
      })
    }
  ) +
  coord_cartesian(clip = "off")+  # 确保文本不会被截断
  # 加粗重复药物名称的标签
  scale_x_discrete(labels = function(drug) {
    sapply(drug, function(x) {
      count <- as.numeric(as.character(top50_drugs_qvalue$Drug_count[top50_drugs_qvalue$Drug == x][1]))
      if (!is.na(count) && count > 1) {
        return(glue::glue("<span style='color:{highlight_color};'><b>{x}</b></span>")) 
      } else {
        return(x)
      }
    })
  }) +
  # 手动添加 qvalue_label 标签到 y 轴上方
  # 手动添加 qvalue_label 标签到 y 轴左侧
  # annotate("text", x = -0.5, y = qvalue_label, label = "0.98", 
  #          fontface = "bold", color = "black", size = 4, hjust = 1)+
  geom_hline(yintercept = qvalue_label, linetype = "dashed", color = "black", size = 1) + 
  annotate("text", x = Inf, y = qvalue_label, label = round(qvalue_label, 2), 
           fontface = "bold", color = "black", size = 4, hjust = -0.1)+

  theme(
    axis.text.x = ggtext::element_markdown(halign = 0.5),  # 解析HTML标签，显示加粗和颜色
    axis.text.y = ggtext::element_markdown(halign = 0)  # 确保 y 轴刻度也能解析加粗标签
  )

# 保存图形
pdf_file <- "Tissue_Drug_GroupedBarplot_Log10_Qvalue_PermP_Stomach_2691.pdf"
ggsave(filename = pdf_file, plot = p, width = 12, height = 10)