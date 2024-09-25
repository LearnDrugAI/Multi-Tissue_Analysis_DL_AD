library(metap)
library(dplyr)
library(data.table)
library(stringr)
library(readxl)

file <- read_excel('KEGG_ATC.xlsx')

DepAnx_drugs <- file %>%
 filter(grepl("^N05B|^N06A", categorySymbols))
# DepAnx_drugs <- file %>%
#  filter(grepl("^N05|^N06", categorySymbols))
# DepAnx_drugs <- file %>%
#   filter(grepl("^N05A", categorySymbols))
# DepAnx_drugs <- file %>%
#   filter(grepl("^N04", categorySymbols))
# # 筛选外部数据集中的Depression_or_Anxiety药物
# DepAnx_drugs <- file %>%
#   filter(grepl("^N03", categorySymbols))

# 假设药物名在第4列，使用 unlist 展开嵌套
DepAnx_drug_list = unique(unlist(DepAnx_drugs[, 4]))##assuming the drug names are listed in the 4th column
# 移除空字符串
DepAnx_drug_list <- DepAnx_drug_list[DepAnx_drug_list != ""]

# 提取随机药物集，并去除空字符串
valid_drug_list <- unique(unlist(file[, 4]))
valid_drug_list <- valid_drug_list[valid_drug_list != ""]
set.seed(123)  # 保证随机性的可重复性
Random_drug_list <- sample(valid_drug_list, size = length(DepAnx_drug_list))
#********************
# 加载重利用数据
#********************
file_names <- list.files(path = "GSCT90042691/transIDtoName", pattern = "^results_5MethodAvg_FiltByP_GCST90042691_.*\\.csv$", full.names = TRUE)
extracted_names <- sub("^GSCT90042691/transIDtoName/results_5MethodAvg_FiltByP_GCST90042691_(.*)\\.csv$", "\\1", file_names)
tissues <- extracted_names
no_tissue = length(tissues)

input_dir <- "GSCT90042691/transIDtoName"
output_dir <- "GSCT90042691/transIDtoName"
# no.drugs=500  # 可调整为你需要的药物数量

pval.oneSamp.t = numeric(no_tissue)
pval.twoSamp.t = numeric(no_tissue)

for (i in 1:no_tissue) {
  CMAP = read.csv(file.path(input_dir, paste0("results_5MethodAvg_FiltByP_GCST90042691_", tissues[i], ".csv")))
  # CMAP <- CMAP[1:no.drugs, ]
  # 去重：基于药物名称去重，保留首次出现的记录
  CMAP <- CMAP %>% distinct(Drug, .keep_all = TRUE)
  # 查找与Depression_or_Anxiety相关的药物在重利用数据中的匹配
  ind.match.trial = grep(paste(DepAnx_drug_list, collapse = "|"), CMAP$Drug, ignore.case = TRUE)
  
  if (length(ind.match.trial) < 1) {
    pval.oneSamp.t[i] <- NA
    pval.twoSamp.t[i] <- NA
    next
  }
  
  # 使用Random_drug_list构建随机药物集
  ind.random.trial = grep(paste(Random_drug_list, collapse = "|"), CMAP$Drug, ignore.case = TRUE)
  
  pval = CMAP$perm.p
  pval[pval == 1] <- 0.999
  pval[pval == 0] <- 1e-4
  
  # 转换为z-scores
  zval = qnorm(pval)
  
  # 与Depression_or_Anxiety药物集合的单样本t检验
  fit1 = t.test(zval[ind.match.trial], alternative = "less", mu = 0, conf.level = 0.95)
  pval.oneSamp.t[i] = fit1$p.value
  
  # 与随机药物集的双样本t检验
  fit2 = t.test(zval[ind.match.trial], zval[ind.random.trial], alternative = "less")
  pval.twoSamp.t[i] = fit2$p.value
}

Fisher_Self <- sumlog(pval.oneSamp.t)
Fisher_Compet <- sumlog(pval.twoSamp.t)
Minp_Self <- minimump(pval.oneSamp.t)
Minp_Compet <- minimump(pval.twoSamp.t)

# 将结果存入数据框
results <- data.frame(
  Meta_SumLog_OneSamp_T = Fisher_Self$p,
  Meta_SumLog_TwoSamp_T = Fisher_Compet$p,
  Meta_MinP_OneSamp_T = Minp_Self$p,
  Meta_MinP_TwoSamp_T = Minp_Compet$p
)

write.csv(results, file.path(output_dir, "GSCT90042691_drug_set_enrichment_results_KEGG_N05BOrN06A_vs_Random.csv"), row.names = FALSE)