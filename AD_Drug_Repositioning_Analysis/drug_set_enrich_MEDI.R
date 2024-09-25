#***********************************
# Drug-set enrichment tests with control for Random Drug Set
#***********************************

library(metap)
library(dplyr)
library(data.table)
library(stringr)

# 读取外部药物数据集
file <- read.csv('MEDI_01212013_HPS.csv', header = TRUE, stringsAsFactors = FALSE, sep = ",")

# 定义疾病列表
Each_Disorder <- c(
  "Anxiety state unspecified",
  "Panic disorder without agoraphobia",
  "Generalized anxiety disorder",
  "Agoraphobia with panic disorder",
  "Agoraphobia without mention of panic attacks",
  "Social phobia",
  "Unspecified nonpsychotic mental disorder",
  "Other and unspecified special symptoms or syndromes; not elsewhere classified",
  "Predominant psychomotor disturbance",
  "Depressive disorder NEC",
  "Nervousness",
  "Psychosexual dysfunction with inhibited sexual excitement",
  "Eating disorder NOS",
  "Persistent disorder of initiating or maintaining sleep",
  "Insomnia; unspecified",
  "Unspecified psychophysiological malfunction",
  "Unspecified hyperkinetic syndrome of childhood",
  "Attention deficit disorder of childhood with hyperactivity",
  "Hyperkinetic syndrome of childhood",
  "Narcolepsy",
  "Restless legs syndrome (RLS)",
  "Alcohol withdrawal delirium",
  "Partial epilepsy without impairment of consciousness without intractable epilepsy",
  "Alcohol dependence syndrome",
  "Paranoid type schizophrenia",
  "Paranoid type schizophrenia; unspecified state",
  "Other and unspecified reactive psychosis",
  "Dementias",
  "Dementia; unspecified; without behavioral disturbance",
  "Major depressive disorder; single episode; unspecified degree",
  "Obsessive-compulsive disorders",
  "Other specified drug-induced mental disorders",
  "Tourette's disorder",
  "Lack of coordination",
  "Other alteration of consciousness",
  "Alcohol withdrawal",
  "Other specified episodic mood disorder",
  "Emotional lability",
  "Catatonic type schizophrenia; unspecified state",
  "Autism disorder",
  "Infantile autism; current or active state",
  "Bulimia nervosa",
  "Opioid abuse",
  "Opioid abuse; unspecified use",
  "Other; mixed; or unspecified drug abuse; unspecified use",
  "Phobia unspecified",
  "Other isolated or specific phobias",
  "Eating disorder NEC",
  "Other specified adjustment reactions",
  "Adjustment reaction with withdrawal",
  "Unspecified sleep disturbance",
  "Unspecified disorder of conduct",
  "Delirium due to conditions classified elsewhere",
  "Other choreas",
  "Convulsions in newborn",
  "Febrile convulsions",
  "Neurosyphilis",
  "Neurosyphilis NOS",
  "Essential and other specified forms of tremor",
  "Altered mental status",
  "Suicidal ideation",
  "Phobic disorders",
  "Undersocialized conduct disorder; unaggressive type; unspecified degree",
  "Schizoaffective disorder",
  "Schizo-affective type schizophrenia; unspecified state",
  "Reading disorder NOS",
  "Unspecified extrapyramidal disease and abnormal movement disorder"
)
Depression_or_Anxiety <- c(
  "Anxiety state unspecified",
  "Panic disorder without agoraphobia",
  "Generalized anxiety disorder",
  "Agoraphobia with panic disorder",
  "Agoraphobia without mention of panic attacks",
  "Social phobia",
  "Depressive disorder NEC",
  "Nervousness",
  "Persistent disorder of initiating or maintaining sleep",
  "Insomnia; unspecified",
  "Major depressive disorder; single episode; unspecified degree",
  "Other specified episodic mood disorder",
  "Emotional lability",
  "Phobia unspecified",
  "Other isolated or specific phobias",
  "Adjustment reaction with withdrawal",
  "Unspecified sleep disturbance",
  "Suicidal ideation",
  "Phobic disorders"
)

All_Psychiatric <- c(
  "Unspecified nonpsychotic mental disorder",
  "Other and unspecified special symptoms or syndromes; not elsewhere classified",
  "Predominant psychomotor disturbance",
  "Psychosexual dysfunction with inhibited sexual excitement",
  "Eating disorder NOS",
  "Unspecified psychophysiological malfunction",
  "Unspecified hyperkinetic syndrome of childhood",
  "Attention deficit disorder of childhood with hyperactivity",
  "Hyperkinetic syndrome of childhood",
  "Narcolepsy",
  "Restless legs syndrome (RLS)",
  "Alcohol withdrawal delirium",
  "Partial epilepsy without impairment of consciousness without intractable epilepsy",
  "Alcohol dependence syndrome",
  "Other and unspecified reactive psychosis",
  "Dementias",
  "Dementia; unspecified; without behavioral disturbance",
  "Obsessive-compulsive disorders",
  "Other specified drug-induced mental disorders",
  "Tourette's disorder",
  "Lack of coordination",
  "Other alteration of consciousness",
  "Alcohol withdrawal",
  "Autism disorder",
  "Infantile autism; current or active state",
  "Bulimia nervosa",
  "Opioid abuse",
  "Opioid abuse; unspecified use",
  "Other; mixed; or unspecified drug abuse; unspecified use",
  "Eating disorder NEC",
  "Other specified adjustment reactions",
  "Unspecified disorder of conduct",
  "Delirium due to conditions classified elsewhere",
  "Other choreas",
  "Convulsions in newborn",
  "Febrile convulsions",
  "Neurosyphilis",
  "Neurosyphilis NOS",
  "Essential and other specified forms of tremor",
  "Altered mental status",
  "Undersocialized conduct disorder; unaggressive type; unspecified degree",
  "Reading disorder NOS",
  "Unspecified extrapyramidal disease and abnormal movement disorder"
)

# 筛选外部数据集中的Depression_or_Anxiety药物
DepAnx_drugs <- file %>%
  filter(INDICATION_DESCRIPTION %in% Depression_or_Anxiety)

DepAnx_drug_list <- unique( as.character(DepAnx_drugs[,2]) )  # 假设药物名在2列中


set.seed(123)  # 保证随机性的可重复性
Random_drug_list <- sample(unique(as.character(file[,2])), size = length(DepAnx_drug_list))

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

write.csv(results, file.path(output_dir, "GSCT90041866_drug_set_enrichment_results_MEDI_All_Psychiatric_vs_Random.csv"), row.names = FALSE)
