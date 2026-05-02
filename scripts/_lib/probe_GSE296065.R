#!/usr/bin/env Rscript
# 轻量探查 GSE296065 Seurat 对象的元数据 / 细胞类型 / 分组
suppressPackageStartupMessages({
  library(Seurat)
})

ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"
RDS  <- file.path(ROOT, "raw_data/00_raw_data/GSE296065_almanzar_flu.rds")

cat(">>> 加载...\n")
t0 <- Sys.time()
seu <- readRDS(RDS)
cat(sprintf(">>> 加载耗时 %.1f 秒\n", as.numeric(Sys.time() - t0, units = "secs")))

cat(sprintf("class: %s\n", paste(class(seu), collapse=",")))
cat(sprintf("ncol (cells): %d, nrow (genes): %d\n", ncol(seu), nrow(seu)))
cat(sprintf("default assay: %s\n", DefaultAssay(seu)))
cat(sprintf("assays: %s\n", paste(Assays(seu), collapse=",")))
if (DefaultAssay(seu) %in% Assays(seu)) {
  cat(sprintf("layers in %s: %s\n",
              DefaultAssay(seu),
              paste(Layers(seu[[DefaultAssay(seu)]]), collapse=",")))
}

cat("\nmetadata 列:\n")
mc <- colnames(seu@meta.data)
for (i in seq_along(mc)) {
  v <- seu@meta.data[[mc[i]]]
  cat(sprintf("  %2d. %-30s  class=%s, n_unique=%d\n",
              i, mc[i], paste(class(v), collapse=","), length(unique(v))))
}

cat("\n=== 候选 cell-type 列的取值 ===\n")
ct_cands <- c("celltype","CellType","cell_type","cell_type_manual",
              "annotation","seurat_clusters","ident","leiden",
              "labels","label","Annotation","cluster","clusters")
for (col in intersect(ct_cands, mc)) {
  cat(sprintf("\n--[%s]--\n", col))
  print(sort(table(seu@meta.data[[col]]), decreasing = TRUE))
}

cat("\n=== 候选 group/treatment 列的取值 ===\n")
grp_cands <- c("group","Group","treatment","Treatment","condition","Condition",
               "sample","Sample","orig.ident","stim","status","disease")
for (col in intersect(grp_cands, mc)) {
  vv <- seu@meta.data[[col]]
  if (length(unique(vv)) <= 30) {
    cat(sprintf("\n--[%s]--\n", col))
    print(sort(table(vv), decreasing = TRUE))
  } else {
    cat(sprintf("\n--[%s]-- (%d unique, 跳过 print)\n", col, length(unique(vv))))
  }
}

# 找包含 Mono/Macroph/Ly6c 关键词的细胞类型 (可能在任意一列)
cat("\n=== 全 metadata 模糊搜 (Mono|Macro|Ly6c) ===\n")
for (col in mc) {
  vv <- seu@meta.data[[col]]
  if (!is.character(vv) && !is.factor(vv)) next
  uv <- unique(as.character(vv))
  hit <- uv[grepl("Mono|Macro|Ly6c", uv, ignore.case = TRUE)]
  if (length(hit) > 0) {
    cat(sprintf("  [%s] -> %s\n", col, paste(hit, collapse=" | ")))
  }
}

# 找 Flu/RTX/Vehicle 关键词
cat("\n=== 全 metadata 模糊搜 (Flu|RTX|Vehicle|Naive|Mock) ===\n")
for (col in mc) {
  vv <- seu@meta.data[[col]]
  if (!is.character(vv) && !is.factor(vv)) next
  uv <- unique(as.character(vv))
  hit <- uv[grepl("Flu|RTX|Vehicle|Naive|Mock|Untreated", uv, ignore.case = TRUE)]
  if (length(hit) > 0) {
    cat(sprintf("  [%s] -> %s\n", col, paste(hit, collapse=" | ")))
  }
}

cat("\n>>> 探查完成\n")
