#!/usr/bin/env Rscript
# 复刻同事 (GSE42639-20260502-1.py) 的方法 + 与正确方法对比
#
# 同事方法 (Python 脚本逻辑):
#   1. 用 detection.tsv (探针 detection p-value 矩阵, 0-1 范围) 作输入
#   2. log2(treat_mean+eps) - log2(control_mean+eps) 当作 logFC
#   3. Welch t-test + BH 校正
#   4. DEG = padj < 0.05 & |log2FC| >= 1
#   5. 没有 normalization, 没有探针 -> 基因, 没有探针级过滤
#
# 我们方法:
#   1. 用 intensity.tsv -> neqc (background + log2 + quantile) -> expr_log2.rds
#   2. 探针级过滤 detection p < 0.05 在 >= 2 样本
#   3. 探针 -> SYMBOL (illuminaMousev2.db, multi-probe -> max-mean rep)
#   4. limma + eBayes
#   5. DEG = adj.P.Val < 0.05 & |logFC| >= 1
#
# 比较:
#   - Mo PR8_100 vs Sham: 同事方法的 DEG vs 我们方法的 DEG, 重叠率
#   - logFC 分布对比
#   - 是否同事方法发现的 "上调" 在我们方法里也是上调

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
})

ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"
source(file.path(ROOT, "scripts/_lib/setup_fonts.R"))

OUT <- file.path(ROOT, "results/celltype_review/colleague_method_check")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# 1. 加载所有 3 种数据
intensity <- read.delim(file.path(ROOT, "processed_data/GSE42639/GSE42639_intensity.tsv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
detection <- read.delim(file.path(ROOT, "processed_data/GSE42639/GSE42639_detection.tsv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.delim(file.path(ROOT, "processed_data/GSE42639/GSE42639_sample_meta.tsv"),
                   stringsAsFactors = FALSE)

# 校正 trim 空格
intensity$probe_id <- trimws(intensity$probe_id)
detection$probe_id <- trimws(detection$probe_id)
stopifnot(identical(intensity$probe_id, detection$probe_id))
sample_cols <- colnames(intensity)[-1]

# Mo Sham vs Mo PR8_100 样本
mo_sham_samples <- meta$sample_id[meta$cell == "Mo" & meta$treat == "Sham"]
mo_pr8_samples  <- meta$sample_id[meta$cell == "Mo" & meta$treat == "100LD50PR8"]
cat(sprintf("Mo Sham n=%d, PR8_100 n=%d\n",
            length(mo_sham_samples), length(mo_pr8_samples)))

# 2. 同事方法: detection p-value t-test
cat("\n>>> 同事方法 (detection p-value t-test)\n")
det_sham <- as.matrix(detection[, mo_sham_samples])
det_pr8  <- as.matrix(detection[, mo_pr8_samples])
det_sham_mean <- rowMeans(det_sham)
det_pr8_mean  <- rowMeans(det_pr8)
det_log2fc <- log2(det_pr8_mean + 1e-10) - log2(det_sham_mean + 1e-10)
det_pval <- vapply(seq_len(nrow(det_sham)), function(i) {
  tryCatch(
    t.test(det_pr8[i,], det_sham[i,], var.equal = FALSE)$p.value,
    error = function(e) 1.0)
}, numeric(1))
det_padj <- p.adjust(det_pval, method = "BH")
colleague_deg_idx <- which(det_padj < 0.05 & abs(det_log2fc) >= 1)
colleague_deg_probes <- detection$probe_id[colleague_deg_idx]
colleague_deg_up   <- detection$probe_id[which(det_padj < 0.05 & det_log2fc >=  1)]
colleague_deg_down <- detection$probe_id[which(det_padj < 0.05 & det_log2fc <= -1)]
cat(sprintf("  同事 DEG: 上 %d / 下 %d / 总 %d (探针级)\n",
            length(colleague_deg_up), length(colleague_deg_down),
            length(colleague_deg_probes)))

# 3. 我们方法: intensity -> neqc -> limma
cat("\n>>> 我们方法 (intensity -> neqc -> limma)\n")
int_M <- as.matrix(intensity[, sample_cols]); rownames(int_M) <- detection$probe_id
det_M <- as.matrix(detection[, sample_cols]); rownames(det_M) <- detection$probe_id
storage.mode(int_M) <- "numeric"; storage.mode(det_M) <- "numeric"

elist <- new("EListRaw", list(E = int_M, other = list(Detection = det_M),
                              targets = meta,
                              genes = data.frame(ProbeID = detection$probe_id)))
elist_n <- neqc(elist)
expr <- elist_n$E
det_filt <- elist_n$other$Detection
keep_probe <- rowSums(det_filt < 0.05) >= 2
expr <- expr[keep_probe, ]

# Mo 子集 + limma
mo_idx <- which(meta$cell == "Mo")
m_mo <- meta[mo_idx, ]; E_mo <- expr[, mo_idx]
treat_R <- c(Sham="Sham","100LD50PR8"="PR8_100","TX91"="TX91",
             "0.6LD50PR8"="PR8_0p6","10LD50PR8"="PR8_10")
m_mo$treat <- factor(treat_R[m_mo$treat],
                     levels=c("Sham","TX91","PR8_0p6","PR8_10","PR8_100"))
design <- model.matrix(~0 + m_mo$treat); colnames(design) <- levels(m_mo$treat)
fit <- lmFit(E_mo, design)
cm <- makeContrasts(PR8_100_vs_Sham = PR8_100 - Sham, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cm), trend = TRUE)
tt <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
tt$probe_id <- rownames(tt)
ours_deg_up   <- tt$probe_id[tt$adj.P.Val < 0.05 & tt$logFC >=  1]
ours_deg_down <- tt$probe_id[tt$adj.P.Val < 0.05 & tt$logFC <= -1]
cat(sprintf("  我们 DEG: 上 %d / 下 %d / 总 %d (探针级)\n",
            length(ours_deg_up), length(ours_deg_down),
            length(ours_deg_up) + length(ours_deg_down)))

# 4. 重叠率
cat("\n>>> 探针级 DEG 重叠分析\n")
overlap_up   <- intersect(colleague_deg_up, ours_deg_up)
overlap_down <- intersect(colleague_deg_down, ours_deg_down)
overlap_any  <- intersect(c(colleague_deg_up, colleague_deg_down),
                          c(ours_deg_up, ours_deg_down))
# 同事 "上调" 在我们里是 "下调" 的 (反向!)
cross_up_to_down <- intersect(colleague_deg_up, ours_deg_down)
cross_down_to_up <- intersect(colleague_deg_down, ours_deg_up)

cat(sprintf("  同事上 ∩ 我们上 (同向): %d / %d (%.1f%%)\n",
            length(overlap_up), length(colleague_deg_up),
            100*length(overlap_up)/max(1,length(colleague_deg_up))))
cat(sprintf("  同事下 ∩ 我们下 (同向): %d / %d (%.1f%%)\n",
            length(overlap_down), length(colleague_deg_down),
            100*length(overlap_down)/max(1,length(colleague_deg_down))))
cat(sprintf("  同事上 ∩ 我们下 (反向!): %d\n", length(cross_up_to_down)))
cat(sprintf("  同事下 ∩ 我们上 (反向!): %d\n", length(cross_down_to_up)))
cat(sprintf("  忽略方向, 任意重叠: %d / 同事总 %d (%.1f%%)\n",
            length(overlap_any), length(c(colleague_deg_up, colleague_deg_down)),
            100*length(overlap_any)/max(1,length(c(colleague_deg_up, colleague_deg_down)))))

# 5. log2FC 分布对比图
cat("\n>>> 生成 logFC 散点对比图\n")
df <- data.frame(
  probe_id   = detection$probe_id,
  colleague_log2FC = det_log2fc,
  colleague_padj   = det_padj
)
ours_lfc_full <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
ours_lfc_full$probe_id <- rownames(ours_lfc_full)
df <- merge(df,
            ours_lfc_full[, c("probe_id","logFC","adj.P.Val")],
            by = "probe_id")
colnames(df)[colnames(df) == "logFC"] <- "ours_log2FC"
colnames(df)[colnames(df) == "adj.P.Val"] <- "ours_padj"

# 4 类标记
df$status <- "ns"
df$status[df$colleague_padj < 0.05 & abs(df$colleague_log2FC) >= 1 &
          df$ours_padj    >= 0.05] <- "仅同事认为是DEG"
df$status[df$ours_padj    < 0.05 & abs(df$ours_log2FC)    >= 1 &
          df$colleague_padj >= 0.05] <- "仅我们认为是DEG"
df$status[df$colleague_padj < 0.05 & abs(df$colleague_log2FC) >= 1 &
          df$ours_padj    < 0.05 & abs(df$ours_log2FC)    >= 1] <- "双方都是DEG"

p <- ggplot(df, aes(colleague_log2FC, ours_log2FC, color = status)) +
  geom_point(size = 0.8, alpha = 0.5) +
  geom_hline(yintercept = c(-1,1), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
  scale_color_manual(values = c("ns" = "grey80",
                                 "仅同事认为是DEG" = "#FF7F0E",
                                 "仅我们认为是DEG" = "#1F77B4",
                                 "双方都是DEG"     = "#D62728")) +
  labs(x = "同事 log2FC (基于 detection p-value)",
       y = "我们 log2FC (基于 neqc 标准化表达)",
       title = "Mo PR8_100 vs Sham: 同事方法 vs 我们方法 logFC 对比",
       subtitle = sprintf("同事 DEG=%d, 我们 DEG=%d, 同向重叠=%d, 反向(矛盾)=%d",
                          length(c(colleague_deg_up,colleague_deg_down)),
                          length(c(ours_deg_up,ours_deg_down)),
                          length(c(overlap_up,overlap_down)),
                          length(c(cross_up_to_down,cross_down_to_up)))) +
  theme_bw(base_size = 11)
ggsave(file.path(OUT, "logFC_compare_scatter.pdf"), p,
       width = 8, height = 7, device = cairo_pdf)
ggsave(file.path(OUT, "logFC_compare_scatter.png"), p,
       width = 8, height = 7, dpi = 200, type = "cairo")

# 6. 写出关键基因诊断: Cx3cr1 探针在两套方法下的位置
cat("\n>>> 关键基因 Cx3cr1 / Klf2 / Arrb1 探针在两套方法的对比\n")
suppressPackageStartupMessages(library(illuminaMousev2.db))
sym_map <- AnnotationDbi::select(illuminaMousev2.db,
                                 keys = detection$probe_id, keytype = "PROBEID",
                                 columns = c("PROBEID","SYMBOL"))
key_genes <- c("Cx3cr1","Klf2","Arrb1","Cd209a","Btla","Il16","Cx3cl1","Gstm1")
for (g in key_genes) {
  pids <- sym_map$PROBEID[sym_map$SYMBOL == g & !is.na(sym_map$SYMBOL)]
  if (length(pids) == 0) next
  for (p in pids) {
    if (!(p %in% detection$probe_id)) next
    coll_lfc <- df$colleague_log2FC[df$probe_id == p]
    coll_padj<- df$colleague_padj  [df$probe_id == p]
    ours_lfc <- df$ours_log2FC     [df$probe_id == p]
    ours_padj<- df$ours_padj       [df$probe_id == p]
    if (length(coll_lfc) == 0 || length(ours_lfc) == 0) next
    cat(sprintf("  %-8s %-15s | 同事 logFC=%+5.2f padj=%.3g | 我们 logFC=%+5.2f padj=%.3g\n",
                g, p, coll_lfc, coll_padj, ours_lfc, ours_padj))
  }
}

# 7. detection p-value 分布特征 (说明为什么 t-test 不合适)
cat("\n>>> detection p 分布特征\n")
det_all <- as.numeric(det_M)
det_lt05 <- sum(det_all < 0.05) / length(det_all)
det_gt95 <- sum(det_all > 0.95) / length(det_all)
cat(sprintf("  全数据点中 detection p<0.05 占比: %.1f%% (强表达探针)\n", 100*det_lt05))
cat(sprintf("  全数据点中 detection p>0.95 占比: %.1f%% (不表达探针)\n", 100*det_gt95))
cat(sprintf("  -> bimodal 分布, t-test 假设的近似高斯分布不成立\n"))

write.table(df, file.path(OUT, "logFC_compare_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(colleague_deg_up = colleague_deg_up,
             colleague_deg_down = colleague_deg_down,
             ours_deg_up = ours_deg_up, ours_deg_down = ours_deg_down,
             overlap = list(up = overlap_up, down = overlap_down,
                            cross_up_to_down = cross_up_to_down,
                            cross_down_to_up = cross_down_to_up),
             df = df),
        file.path(OUT, "colleague_vs_ours.rds"))

cat(sprintf("\n>>> 输出: %s\n", OUT))
