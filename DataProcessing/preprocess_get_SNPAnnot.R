library(VariantAnnotation)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(GenomicRanges)

# 文件路径和存在性检查
base_file <- "NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_"
vcf_suffix <- ".recalibrated_variants.Broad_Rush.vcf.gz"

# 设置 VCF 文件读取参数和分块大小
param <- ScanVcfParam(fixed = c("ALT"), info = c("AC"), geno = NA)  # 使用AC字段
chunksize <- 100000  # 将块大小调整为100,000

# 本地注释数据
genome <- BSgenome.Hsapiens.UCSC.hg19
all_snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37

# 循环处理每个染色体
for (chr in 1:22) {
  file <- paste0(base_file, chr, vcf_suffix)
  vcf.path <- file.path(file)
  if (!file.exists(vcf.path)) {
    cat("VCF file does not exist:", vcf.path, "\n")
    next
  }
  
  # 打开 VCF 文件
  myVcfFile <- TabixFile(vcf.path, index = paste0(vcf.path, ".tbi"), yieldSize = chunksize)
  open(myVcfFile)
  
  # 初始化注释文件
  annotation_file <- paste0("snp_annot.chr", chr, ".txt")
  header_written <- FALSE
  
  # 处理每个块
  while (TRUE) {
    vcf_chunk <- readVcf(myVcfFile, param = param)
    if (nrow(vcf_chunk) == 0) break
    
    # 过滤非SNVs和多个ALT等位基因
    snvs <- isSNV(vcf_chunk, singleAltOnly = TRUE)
    vcf_chunk <- vcf_chunk[snvs]
    
    # 提取所需信息
    snp_info <- as.data.table(as.data.frame(rowRanges(vcf_chunk)))
    snp_info <- snp_info[, .(seqnames, start, REF, ALT)]
    setnames(snp_info, c("seqnames", "start", "REF", "ALT"), c("chromosome", "Pos", "ref_vcf", "alt_vcf"))
    
    # 确保 ALT 列为字符向量
    snp_info[, alt_vcf := sapply(alt_vcf, function(x) as.character(x[[1]]))]
    
    # 添加 snp_id_originalVCF 和 VariantID 列
    snp_info[, snp_id_originalVCF := paste0("snp_", chromosome, "_", Pos)]
    snp_info[, varID := paste0(chromosome, ":", Pos, "_", ref_vcf, "/", alt_vcf)]
    snp_info[, Num_alt_per_site := unlist(info(vcf_chunk)$AC)]
    
    # 构建 GRanges 对象并进行本地注释
    positions <- GRanges(seqnames = snp_info$chromosome, ranges = IRanges(start = snp_info$Pos, width = 1))
    
    # 本地注释
    my_snps <- snpsByOverlaps(all_snps, positions)
    
    # 将注释信息添加到 snp_info 中
    rsid_map <- data.table(
      chr_pos = paste0(seqnames(my_snps), ":", start(my_snps)),
      rsid = my_snps$RefSNP_id
    )
    snp_info[, rsid := rsid_map$rsid[match(paste0(chromosome, ":", Pos), rsid_map$chr_pos)]]
    
    # 选择并排序所需列
    snp_info <- snp_info[, .(chromosome, Pos, varID, ref_vcf, alt_vcf, snp_id_originalVCF, rsid, Num_alt_per_site)]
    
    # 去除 rsid 列为空的行
    snp_annot <- snp_info[!is.na(rsid) & rsid != ""]
    
    # 写入注释文件
    if (!header_written) {
      fwrite(snp_annot, file = annotation_file, sep = "\t", col.names = TRUE, row.names = FALSE)
      header_written = TRUE
    } else {
      fwrite(snp_annot, file = annotation_file, sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
    }
    
    cat("Processed lines:", nrow(vcf_chunk), "\n")
    cat("Annotation file written:", file.exists(annotation_file), "\n")
    
    # 强制执行垃圾回收，清理内存
    rm(vcf_chunk, snvs, snp_info, positions, my_snps, rsid_map)
    gc()
  }
  
  # 关闭 VCF 文件
  close(myVcfFile)
}

cat("All annotation files generated.\n")
