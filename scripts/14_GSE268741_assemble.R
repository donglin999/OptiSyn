#!/usr/bin/env Rscript
# GSE268741 (Gannot et al 2024 NTS scRNA-seq) 数据组装
# 将 31 个 xlsx triplet 切片合并 -> sparse Matrix -> Seurat 对象 (cache)
#
# 输入: raw_data/00_raw_data/GSE268741/
#   GSM8298251_barcodes.xlsx     (9638 cells)
#   GSM8298251_features.xlsx     (33696 genes; ENSEMBL + SYMBOL)
#   GSM8298251_matrix_part01-31.xlsx  (32,201,943 triplets)
#   GSM8298251_matrix_meta.xlsx
#
# 输出: processed_data/GSE268741/
#   counts_sparse.rds    dgCMatrix (gene SYMBOL x cell barcode)
#   features.tsv         feature_id / feature_name / feature_type
#   barcodes.tsv

suppressPackageStartupMessages({
  library(readxl)
  library(Matrix)
})

ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"
RDIR <- file.path(ROOT, "raw_data/00_raw_data/GSE268741")
PROC <- file.path(ROOT, "processed_data/GSE268741")
dir.create(PROC, recursive = TRUE, showWarnings = FALSE)

OUT_RDS <- file.path(PROC, "counts_sparse.rds")
if (file.exists(OUT_RDS)) {
  cat(sprintf(">>> %s 已存在, 跳过组装\n", OUT_RDS))
  quit(status = 0)
}

t0 <- Sys.time()

# 1. barcodes + features
cat(">>> 读 barcodes / features\n")
bc <- as.data.frame(read_excel(file.path(RDIR, "GSM8298251_barcodes.xlsx")))
ft <- as.data.frame(read_excel(file.path(RDIR, "GSM8298251_features.xlsx")))
cat(sprintf("    barcodes: %d, features: %d\n", nrow(bc), nrow(ft)))

# 2. 31 个 triplet 切片
mat_files <- sort(list.files(RDIR, pattern = "^GSM8298251_matrix_part[0-9]+\\.xlsx$",
                              full.names = TRUE))
cat(sprintf(">>> 读 %d 个 triplet xlsx 切片\n", length(mat_files)))

trips <- vector("list", length(mat_files))
for (i in seq_along(mat_files)) {
  ti <- Sys.time()
  d <- read_excel(mat_files[i])
  trips[[i]] <- d
  cat(sprintf("    [%2d/%d] %s rows=%d (%.1fs)\n",
              i, length(mat_files),
              basename(mat_files[i]), nrow(d),
              as.numeric(Sys.time() - ti, units = "secs")))
}
trip <- do.call(rbind, trips)
rm(trips); gc()
cat(sprintf(">>> 合并 triplet: %d rows\n", nrow(trip)))

# 3. 构造稀疏矩阵 (gene x cell)
i_idx <- as.integer(trip$feature_index_1based)
j_idx <- as.integer(trip$barcode_index_1based)
v_val <- as.numeric(trip$count)
M <- sparseMatrix(i = i_idx, j = j_idx, x = v_val,
                  dims = c(nrow(ft), nrow(bc)))
rownames(M) <- make.unique(ft$feature_name)
colnames(M) <- bc$barcode
cat(sprintf(">>> 稀疏矩阵: %d gene x %d cell, nnz=%d\n",
            nrow(M), ncol(M), length(M@x)))

# 4. 保存
saveRDS(M, OUT_RDS)
write.table(ft, file.path(PROC, "features.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bc, file.path(PROC, "barcodes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf(">>> 组装耗时: %.1f 分钟\n",
            as.numeric(Sys.time() - t0, units = "mins")))
cat(sprintf(">>> 输出: %s\n", PROC))
