#!/usr/bin/env Rscript
# GSE31022 预处理: neqc 标准化 + 探针 -> 基因
#
# 数据: raw_data/00_raw_data/GSE31022_non-normalized.txt
#   首列: PROBE_ID (Illumina ILMN_*)
#   其后 21 个样本, 每样本 2 列 (intensity, Detection Pval)
#   样本: Ctrl1/2/3, D{1..6}-{1..3}  (3 Ctrl + 6 day x 3 rep = 21)
#
# 输出:
#   processed_data/GSE31022/expr_log2.rds        (gene x sample)
#   processed_data/GSE31022/sample_meta.rds      (data.frame)
#   processed_data/GSE31022/probe2gene.tsv

suppressPackageStartupMessages({
  library(limma)
  library(illuminaMousev2.db)
})

ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"
RAW  <- file.path(ROOT, "raw_data/00_raw_data/GSE31022_non-normalized.txt")
OUT  <- file.path(ROOT, "processed_data/GSE31022")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cat(">>> 读取原始矩阵\n")
df <- read.delim(RAW, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(colnames(df)[1] == "PROBE_ID")
ncol_total <- ncol(df)
stopifnot((ncol_total - 1) %% 2 == 0)

# 拆 intensity / detection (交错列)
n_pairs <- (ncol_total - 1) / 2
int_idx <- 1 + 2 * (seq_len(n_pairs) - 1) + 1   # 2,4,6,... -> intensity
det_idx <- 1 + 2 * (seq_len(n_pairs) - 1) + 2   # 3,5,7,... -> detection
sample_names <- colnames(df)[int_idx]
cat(sprintf("    样本数: %d\n", length(sample_names)))
cat("    样本: ", paste(sample_names, collapse=", "), "\n")

probes <- trimws(df$PROBE_ID)
E_raw  <- as.matrix(df[, int_idx, drop = FALSE]); colnames(E_raw) <- sample_names
P_det  <- as.matrix(df[, det_idx, drop = FALSE]); colnames(P_det) <- sample_names
rownames(E_raw) <- probes; rownames(P_det) <- probes
storage.mode(E_raw) <- "numeric"
storage.mode(P_det) <- "numeric"
cat(sprintf("    raw matrix: %d 探针 x %d 样本\n", nrow(E_raw), ncol(E_raw)))

# 构造样本 meta
parse_sample <- function(s) {
  if (grepl("^Ctrl", s)) {
    list(group = "Ctrl", day = 0L,  rep = as.integer(sub("Ctrl", "", s)))
  } else {
    mt <- regmatches(s, regexec("^D([0-9]+)-([0-9]+)$", s))[[1]]
    list(group = paste0("D", mt[2]), day = as.integer(mt[2]), rep = as.integer(mt[3]))
  }
}
parsed <- lapply(sample_names, parse_sample)
meta <- data.frame(
  sample_id = sample_names,
  group     = vapply(parsed, function(x) as.character(x$group), character(1)),
  day       = vapply(parsed, function(x) as.integer(x$day),    integer(1)),
  rep       = vapply(parsed, function(x) as.integer(x$rep),    integer(1)),
  stringsAsFactors = FALSE
)
rownames(meta) <- meta$sample_id
cat(">>> 样本分布:\n"); print(table(meta$group))

# neqc 标准化
elist <- new("EListRaw",
             list(E = E_raw, other = list(Detection = P_det),
                  targets = meta, genes = data.frame(ProbeID = probes)))
cat(">>> neqc 标准化\n")
elist_n <- neqc(elist)
expr <- elist_n$E
det  <- elist_n$other$Detection
cat(sprintf("    标准化后矩阵: %d x %d\n", nrow(expr), ncol(expr)))

# 探针表达过滤: detection p < 0.05 在至少 2 个样本
keep_probe <- rowSums(det < 0.05) >= 2
cat(sprintf(">>> 探针检测过滤 detP<0.05 >= 2 样本: 保留 %d / %d\n",
            sum(keep_probe), length(keep_probe)))
expr <- expr[keep_probe, , drop = FALSE]
det  <- det[keep_probe, , drop = FALSE]

# 探针 -> 基因 (illuminaMousev2.db)
cat(">>> 探针 -> 基因符号映射\n")
probe_keep <- rownames(expr)
sym_map <- AnnotationDbi::select(illuminaMousev2.db,
                                 keys = probe_keep, keytype = "PROBEID",
                                 columns = c("PROBEID","SYMBOL","ENTREZID"))
sym_map <- sym_map[!is.na(sym_map$SYMBOL), ]
# 一对多探针: 取第一个 SYMBOL
sym_map <- sym_map[!duplicated(sym_map$PROBEID), ]
write.table(sym_map, file.path(OUT, "probe2gene.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("    带 SYMBOL 的探针: %d\n", nrow(sym_map)))

# 多探针映射同基因: 取每基因表达均值最高的探针 (经典做法)
expr2 <- expr[sym_map$PROBEID, , drop = FALSE]
gene_v <- sym_map$SYMBOL
ave_per_probe <- rowMeans(expr2)
ord <- order(gene_v, -ave_per_probe)
expr_ord <- expr2[ord, , drop = FALSE]
gene_ord <- gene_v[ord]
keep_first <- !duplicated(gene_ord)
expr_gene <- expr_ord[keep_first, , drop = FALSE]
rownames(expr_gene) <- gene_ord[keep_first]
cat(sprintf(">>> 基因级矩阵: %d 基因 x %d 样本\n",
            nrow(expr_gene), ncol(expr_gene)))

# 保存
saveRDS(expr_gene, file.path(OUT, "expr_log2.rds"))
saveRDS(meta,      file.path(OUT, "sample_meta.rds"))
cat(sprintf(">>> 输出: %s\n", OUT))
