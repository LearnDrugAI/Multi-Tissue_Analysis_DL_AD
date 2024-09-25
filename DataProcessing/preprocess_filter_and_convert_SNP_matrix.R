library(VariantAnnotation)
library(dplyr)
library(data.table)

# 文件路径和存在性检查
file <- "NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_3.recalibrated_variants.Broad_Rush.vcf.gz"
vcf.path <- file.path(file)
if (!file.exists(vcf.path)) {
  stop("VCF file does not exist.")
}

# 读取样本标识符
samples.vcf <- samples(scanVcfHeader(vcf.path))
length(samples.vcf)  # 1196

# 读取 biospecimen 和 clinical 数据
tbl.biospecimen <- read.table("ROSMAP_biospecimen_metadata.csv", sep=",", header=TRUE, as.is=TRUE)
tbl.clinical <- read.table("ROSMAP_clinical.csv", sep=",", as.is=TRUE, header=TRUE, nrow=-1)

# 合并样本标识符和 biospecimen 数据
tbl.ids <- data.frame(vcf=samples.vcf, stringsAsFactors=FALSE)
tbl.ids2 <- merge(tbl.ids, tbl.biospecimen[, c("specimenID", "organ", "tissue", "individualID")], by.x="vcf", by.y="specimenID", all.x=TRUE)
tbl.ids3 <- merge(tbl.ids2, tbl.clinical, by="individualID", all.x=TRUE)
dim(tbl.ids3)

# 读取 RNA 表达数据
mtx.rna <- get(load("rosmap.14235x632.RData"))
samples.rna <- colnames(mtx.rna)
length(samples.rna)  # 632

# 合并 VCF 和 RNA 样本标识符
samples.vcf.rna <- c(samples.vcf, samples.rna)
samples.all <- intersect(samples.vcf.rna, tbl.biospecimen$specimenID)
selected_samples <- samples.vcf[samples.vcf %in% samples.all]

# 设置 VCF 文件读取参数和分块大小
param <- ScanVcfParam(fixed = c("ALT"), geno = "GT", info = NA)
chunksize <- 500000

# 打开 VCF 文件
myVcfFile <- TabixFile(vcf.path, index = paste0(vcf.path, ".tbi"), yieldSize = chunksize)
open(myVcfFile)

# 初始化文件并写入头部
csvFileName <- sub(".vcf.gz$", ".snptable.csv", file)
header_written <- FALSE

while (TRUE) {
  vcf_chunk <- readVcf(myVcfFile, param = param)
  if (nrow(vcf_chunk) == 0) break
  
  # 过滤非SNVs和多个ALT等位基因
  snvs <- isSNV(vcf_chunk, singleAltOnly = TRUE)
  vcf_chunk <- vcf_chunk[snvs]
  
  # 转换为SNP矩阵
  snps <- genotypeToSnpMatrix(vcf_chunk)
  z <- apply(snps$genotypes@.Data, 2, function(x) as.integer(x))
  rownames(z) <- colnames(geno(vcf_chunk)$GT)
  snps_df <- as.data.frame(z)
  selected_snps <- snps_df[selected_samples, , drop = FALSE] 
  
  # 写入CSV文件
  if (!header_written) {
    fwrite(selected_snps, file = csvFileName, sep = ",", col.names = TRUE, row.names = TRUE)
    header_written <- TRUE
  } else {
    fwrite(selected_snps, file = csvFileName, sep = ",", col.names = FALSE, row.names = TRUE, append = TRUE)
  }
  
  cat("Processed lines:", nrow(vcf_chunk), "\n")
}

# 关闭 VCF 文件
close(myVcfFile)
