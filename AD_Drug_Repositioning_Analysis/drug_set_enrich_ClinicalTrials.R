#***********************************
# Drug-set enrichment tests with control for Random Drug Set
#***********************************

library(metap)
library(dplyr)
library(data.table)
library(stringr)

file <- read.csv('DrugBank-DO.tsv', header = TRUE, stringsAsFactors = FALSE, sep = "\t")
Irritability <- c(
  "anxiety disorder",
  "autism disorder",
  "attention deficit hyperactivity disorder",
  "bipolar disorder", 
  "depression", 
  "intermittent explosive disorder",
  "borderline personality disorder",
  "obsessive-compulsive disorder", 
  "post-traumatic stress disorder",
  "schizophrenia",
  "substance abuse",
  "sleep disorders",
  "thyroid disorders",
  "dementia",
  "autism disorder",
  "attention deficit hyperactivity disorder",
  "bipolar disorder",
  "depression",
  "intermittent explosive disorder",
  "borderline personality disorder",
  "obsessive-compulsive disorder",
  "post-traumatic stress disorder",
  "schizophrenia",
  "substance abuse",
  "thyroid disorders",
  "dementia",
  "multiple sclerosis",
  "restless legs syndrome",
  "amyotrophic lateral sclerosis",
  "endometriosis",
  "autism disorder",
  "attention deficit hyperactivity disorder",
  "bipolar disorder",
  "depression",
  "intermittent explosive disorder",
  "borderline personality disorder",
  "obsessive-compulsive disorder",
  "post-traumatic stress disorder",
  "schizophrenia",
  "substance abuse",
  "thyroid disorders",
  "dementia",
  "multiple sclerosis",
  "restless legs syndrome",
  "amyotrophic lateral sclerosis",
  "endometriosis",
  "schizoaffective disorder",
  "somatoform disorder",
  "fibromyalgia",
  "prion disease",
  "postpartum depression",
  "chronic fatigue syndrome",
  "chronic myeloid leukemia",
  "myasthenia gravis",
  "eating disorder",
  "mycosis fungoides",
  "anorexia nervosa"
)

Worrier <- c(
  "anxiety disorder",
  "generalized anxiety disorder",
  "panic disorder",
  "social phobia",
  "obsessive-compulsive disorder",
  "post-traumatic stress disorder",
  "depression",
  "specific phobias",
  "schizophrenia",
  "substance abuse",
  "eating disorders",
  "thyroid disorders",
  "irritable bowel syndrome",
  "anxiety disorder",
  "generalized anxiety disorder",
  "panic disorder",
  "social phobia",
  "obsessive-compulsive disorder",
  "post-traumatic stress disorder",
  "depression",
  "specific phobias",
  "schizophrenia",
  "substance abuse",
  "eating disorders",
  "thyroid disorders",
  "irritable bowel syndrome",
  "panic disorder",
  "panic disorder",
  "cystic fibrosis",
  "anxiety disorder",
  "generalized anxiety disorder",
  "panic disorder",
  "social phobia",
  "obsessive-compulsive disorder",
  "post-traumatic stress disorder",
  "depression",
  "specific phobias",
  "schizophrenia",
  "substance abuse",
  "eating disorders",
  "thyroid disorders",
  "irritable bowel syndrome",
  "cystic fibrosis",
  "schizoaffective disorder",
  "panic disorder",
  "agoraphobia",
  "hypochondriasis",
  "somatoform disorder",
  "phobic disorder",
  "neurotic disorder",
  "insomnia",
  "somatization disorder",
  "panic disorder",
  "Gilles de la Tourette syndrome",
  "dysthymic disorder",
  "dissociative disorder",
  "borderline personality disorder",
  "impulse control disorder"
)

Alzheimer <- c(
  "Alzheimer's disease",
  "dementia",
  "Lewy body dementia",
  "vascular dementia",
  "mild cognitive impairment",
  "Parkinson's disease",
  "Huntington's disease",
  "amyotrophic lateral sclerosis",
  "frontotemporal dementia",
  "Creutzfeldt-Jakob disease",
  "neurodegenerative diseases",
  "traumatic brain injury",
  "Down syndrome",
  "Alzheimer's disease",
  "dementia",
  "Lewy body dementia",
  "vascular dementia",
  "mild cognitive impairment",
  "Parkinson's disease",
  "Huntington's disease",
  "amyotrophic lateral sclerosis",
  "frontotemporal dementia",
  "Creutzfeldt-Jakob disease",
  "neurodegenerative diseases",
  "traumatic brain injury",
  "Down syndrome",
  "Pick's disease",
  "cognitive disorder",
  "dementia",
  "amyotrophic lateral sclerosis",
  "Alzheimer's disease",
  "dementia",
  "Lewy body dementia",
  "vascular dementia",
  "mild cognitive impairment",
  "Parkinson's disease",
  "Huntington's disease",
  "amyotrophic lateral sclerosis",
  "frontotemporal dementia",
  "Creutzfeldt-Jakob disease",
  "neurodegenerative diseases",
  "traumatic brain injury",
  "Down syndrome",
  "Pick's disease",
  "cognitive disorder",
  "amyotrophic lateral sclerosis",
  "primary cerebellar degeneration",
  "neurofibromatosis",
  "neuroblastoma",
  "Wernicke-Korsakoff syndrome",
  "progressive supranuclear palsy"
)


# 筛选外部数据集中的Depression_or_Anxiety药物
DepAnx_drugs <- file %>%
  filter(doid_name %in% Alzheimer)
DepAnx_drug_list <-unique( as.character(DepAnx_drugs[,3]) )  # 假设药物名在"DRUG_NAME"列中
# 移除空字符串
DepAnx_drug_list <- DepAnx_drug_list[DepAnx_drug_list != ""]
# 转义特殊字符
escape_special_chars <- function(drug) {
  str_replace_all(drug, "([.|(){}\\[\\]^$+*?\\\\])", "\\\\\\1")
}

escaped_drugs <- sapply(DepAnx_drug_list, escape_special_chars)

# 创建正则表达式
pattern1 <- paste(escaped_drugs, collapse = "|")
# 提取随机药物集，并去除空字符串
valid_drug_list <- unique(as.character(file[, 3]))
valid_drug_list <- valid_drug_list[valid_drug_list != ""]

# 从外部数据集中随机抽取与DepAnx_drug_list相同数量的药物
set.seed(123)  # 保证随机性的可重复性
Random_drug_list <- sample(valid_drug_list, size = length(DepAnx_drug_list))

# 转义随机药物集中的特殊字符
escaped_drugs_random <- sapply(Random_drug_list, escape_special_chars)
pattern2 <- paste(escaped_drugs_random, collapse = "|")

#********************
# 加载重利用数据
#********************
file_names <- list.files(path = "GSCT90041866/transIDtoName", pattern = "^results_5MethodAvg_FiltByP_GCST90041866_.*\\.csv$", full.names = TRUE)
extracted_names <- sub("^GSCT90041866/transIDtoName/results_5MethodAvg_FiltByP_GCST90041866_(.*)\\.csv$", "\\1", file_names)
tissues <- extracted_names
no_tissue = length(tissues)

input_dir <- "GSCT90041866/transIDtoName"
output_dir <- "GSCT90041866/transIDtoName"
# no.drugs=500  # 可调整为你需要的药物数量

pval.oneSamp.t = numeric(no_tissue)
pval.twoSamp.t = numeric(no_tissue)

for (i in 1:no_tissue) {
  CMAP = read.csv(file.path(input_dir, paste0("results_5MethodAvg_FiltByP_GCST90041866_", tissues[i], ".csv")))
  # CMAP <- CMAP[1:no.drugs, ]
  # 去重：基于药物名称去重，保留首次出现的记录
  CMAP <- CMAP %>% distinct(Drug, .keep_all = TRUE)
  # 查找与Depression_or_Anxiety相关的药物在重利用数据中的匹配
  ind.match.trial = grep( pattern1, CMAP$Drug  , ignore.case=TRUE)
  if (length(ind.match.trial) < 1) {
    pval.oneSamp.t[i] <- NA
    pval.twoSamp.t[i] <- NA
    next
  }
  
  # 使用Random_drug_list构建随机药物集
  ind.random.trial = grep( pattern2, CMAP$Drug  , ignore.case=TRUE)
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

write.csv(results, file.path(output_dir, "GSCT90041866_drug_set_enrichment_results_ClinicalTrials_Alzheimer_vs_Random.csv"), row.names = FALSE)
