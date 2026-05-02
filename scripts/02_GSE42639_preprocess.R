#!/usr/bin/env Rscript
# GSE42639 预处理：limma::neqc 标准化 + 探针→基因 + 检测 p 过滤
# 输入: processed_data/GSE42639/{intensity,detection,sample_meta}.tsv
# 输出: processed_data/GSE42639/expr_log2.rds        (gene × sample, log2 quantile-normalized)
#       processed_data/GSE42639/sample_meta.rds      (data.frame，列与 expr 同序)
#       processed_data/GSE42639/probe2gene.tsv       (探针注释参考)

suppressPackageStartupMessages({
  library(limma)
  library(illuminaMousev2.db)
})

ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"
INDIR <- file.path(ROOT, "processed_data/GSE42639")

cat(">>> 读取 intensity / detection / meta\n")
int_df <- read.delim(file.path(INDIR, "GSE42639_intensity.tsv"),
                     check.names = FALSE, stringsAsFactors = FALSE)
det_df <- read.delim(file.path(INDIR, "GSE42639_detection.tsv"),
                     check.names = FALSE, stringsAsFactors = FALSE)
meta   <- read.delim(file.path(INDIR, "GSE42639_sample_meta.tsv"),
                     check.names = FALSE, stringsAsFactors = FALSE)

# 探针 ID 末尾有空格，统一去除
int_df$probe_id <- trimws(int_df$probe_id)
det_df$probe_id <- trimws(det_df$probe_id)
stopifnot(identical(int_df$probe_id, det_df$probe_id))
stopifnot(identical(colnames(int_df), colnames(det_df)))
probes <- int_df$probe_id

# 列顺序与 meta 对齐
sample_cols <- colnames(int_df)[-1]
stopifnot(setequal(sample_cols, meta$sample_id))
meta <- meta[match(sample_cols, meta$sample_id), , drop = FALSE]
rownames(meta) <- meta$sample_id

E_raw <- as.matrix(int_df[, -1]);  rownames(E_raw) <- probes
P_det <- as.matrix(det_df[, -1]);  rownames(P_det) <- probes
storage.mode(E_raw) <- "numeric"
storage.mode(P_det) <- "numeric"
cat(sprintf("    raw matrix: %d 探针 × %d 样本\n", nrow(E_raw), ncol(E_raw)))

# 构造 EListRaw 给 limma::neqc(Illumina BeadChip + Detection p)
# neqc: background correction (normexp) + log2 + quantile normalization
elist <- new("EListRaw", list(E = E_raw, other = list(Detection = P_det),
                              targets = meta, genes = data.frame(ProbeID = probes)))
cat(">>> neqc 标准化（normexp 背景校正 + log2 + quantile）\n")
elist_n <- neqc(elist)
expr <- elist_n$E
det  <- elist_n$other$Detection
cat(sprintf("    标准化后矩阵: %d × %d\n", nrow(expr), ncol(expr)))

# 探针级表达过滤：每个处理组（cell × treat）最小重复数为 2，
# 因此要求至少 2 个样本 detection p < 0.05 才保留
keep_probe <- rowSums(det < 0.05) >= 2
cat(sprintf(">>> 探针检测过滤 detP<0.05 ≥ 2 样本: 保留 %d / %d\n",
            sum(keep_probe), length(keep_probe)))
expr <- expr[keep_probe, , drop = FALSE]
det  <- det[keep_probe, , drop = FALSE]

# 探针 → 基因符号 注释
cat(">>> illuminaMousev2.db 注释\n")
ann <- AnnotationDbi::select(illuminaMousev2.db,
                             keys = rownames(expr),
                             keytype = "PROBEID",
                             columns = c("SYMBOL", "ENTREZID")) |>
  unique()
ann <- ann[!is.na(ann$SYMBOL) & ann$SYMBOL != "", , drop = FALSE]
# 若一探针映射多个 SYMBOL，仅保留第一个
ann <- ann[!duplicated(ann$PROBEID), , drop = FALSE]
write.table(ann, file.path(INDIR, "probe2gene.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

expr2 <- expr[rownames(expr) %in% ann$PROBEID, , drop = FALSE]
ann2  <- ann[match(rownames(expr2), ann$PROBEID), , drop = FALSE]
cat(sprintf("    可注释探针: %d\n", nrow(expr2)))

# 同基因多探针：保留平均表达最高的探针
avg <- rowMeans(expr2)
ord <- order(ann2$SYMBOL, -avg)
expr2 <- expr2[ord, , drop = FALSE]
ann2  <- ann2[ord, , drop = FALSE]
keep_top <- !duplicated(ann2$SYMBOL)
expr_gene <- expr2[keep_top, , drop = FALSE]
rownames(expr_gene) <- ann2$SYMBOL[keep_top]
cat(sprintf(">>> 基因级矩阵: %d 基因 × %d 样本\n",
            nrow(expr_gene), ncol(expr_gene)))

# 保存
saveRDS(expr_gene, file.path(INDIR, "expr_log2.rds"))
saveRDS(meta,      file.path(INDIR, "sample_meta.rds"))
cat(">>> 保存:\n  expr_log2.rds\n  sample_meta.rds\n  probe2gene.tsv\n")

# 简单 QC：每细胞类型每处理组的样本数
cat("\n=== sample_meta cell × treat ===\n")
print(addmargins(table(meta$cell, meta$treat)))
cat("\n=== expr range ===\n")
print(summary(as.vector(expr_gene[sample(nrow(expr_gene), 5000), ])))
