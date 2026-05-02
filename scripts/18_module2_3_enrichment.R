#!/usr/bin/env Rscript
# 模块 2/3: Tac1+ 神经元流感信号富集 + 跨数据集通路一致性
#
# 模块 2 流程:
#   2.1 Tac1+Vglut2+ vs 其它 NTS 神经元 FindMarkers (定义 Tac1+ 特异 marker)
#   2.2 Tac1+ 特异 marker 的 GO 富集 (limma::goana, biological process)
#   2.3 用 GSE161878 上调 DEG 144 在 Tac1+ 中做 enrichment 检验 (hypergeometric)
#   2.4 核心靶点 (Cx3cr1, Gpr151, Gstm1, Cx3cl1) 在 Tac1+ 共定位 VlnPlot
#
# 模块 3 流程:
#   3.1 GSE161878 上调 DEG 的 GO 富集
#   3.2 Tac1+ 特异 marker 的 GO 富集
#   3.3 共享通路韦恩图 + 热图 (按 Top 通路对齐)
#   3.4 核心信号传导分子清单 (vagus 上调 ∩ Tac1+ marker)
#
# 输出: results/GSE268741/module2_3_enrichment/

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(limma)
  library(org.Mm.eg.db)
})

set.seed(42)
TS <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"

# 启用中文字体 (showtext + macOS STHeiti)
source(file.path(ROOT, "scripts/_lib/setup_fonts.R"))
OUT  <- file.path(ROOT, "results/GSE268741/module2_3_enrichment")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# 1. 加载 ------------------------------------------------------------------
cat(">>> 加载\n")
seu <- readRDS(file.path(ROOT, "processed_data/GSE268741/seurat_module1.rds"))
strict <- read.delim(file.path(ROOT,
            "results/Mo_protective/strict/core_targets_strict.tsv"),
            stringsAsFactors = FALSE)
genes58 <- strict$gene_symbol

deg_full <- read.delim(file.path(ROOT, "results/GSE161878/GSE161878_DEA_full.tsv"),
                       stringsAsFactors = FALSE)
deg_full <- deg_full[!is.na(deg_full$SYMBOL), ]
flu_up   <- deg_full$SYMBOL[!is.na(deg_full$FDR) & deg_full$FDR < 0.05 &
                            deg_full$log2FC_shrink >= 0.585]
flu_down <- deg_full$SYMBOL[!is.na(deg_full$FDR) & deg_full$FDR < 0.05 &
                            deg_full$log2FC_shrink <= -0.585]
cat(sprintf("    GSE161878 DEG (剔 D-suffix 主分析): up=%d, down=%d\n",
            length(flu_up), length(flu_down)))

# 2. 定义 Tac1+Vglut2+ vs 其它神经元 ---------------------------------------
expr_data <- GetAssayData(seu, assay = "RNA", layer = "data")
neuron_panel <- intersect(c("Snap25","Stmn2","Syp","Map2","Syn1","Tubb3",
                            "Rbfox3","Eno2","Uchl1"), rownames(expr_data))
neuron_score <- Matrix::colMeans(expr_data[neuron_panel, ] > 0)
tac1_pos    <- expr_data["Tac1", ] > 0
slc17a6_pos <- expr_data["Slc17a6", ] > 0

is_neuron <- neuron_score >= 0.5
is_tac1   <- is_neuron & tac1_pos & slc17a6_pos
is_other_n<- is_neuron & !is_tac1
seu$nlabel <- "non_neuron"
seu$nlabel[is_tac1]   <- "Tac1+Vglut2+"
seu$nlabel[is_other_n] <- "Other_neuron"
cat("    分组:\n"); print(table(seu$nlabel))

# 限制对比在神经元内
seu_neu <- subset(seu, cells = colnames(seu)[seu$nlabel %in% c("Tac1+Vglut2+","Other_neuron")])
Idents(seu_neu) <- "nlabel"

# 3. 模块 2.1 FindMarkers -------------------------------------------------
cat(">>> 模块 2.1 Tac1+Vglut2+ vs Other_neuron FindMarkers\n")
mk <- FindMarkers(seu_neu,
                  ident.1 = "Tac1+Vglut2+",
                  ident.2 = "Other_neuron",
                  test.use = "wilcox",
                  logfc.threshold = 0.25,
                  min.pct = 0.10,
                  only.pos = FALSE,
                  verbose = FALSE)
mk$gene_symbol <- rownames(mk); rownames(mk) <- NULL
mk <- mk[, c("gene_symbol","avg_log2FC","p_val","p_val_adj","pct.1","pct.2")]
colnames(mk) <- c("gene_symbol","avg_log2FC","p_val","p_val_adj",
                  "pct_Tac1","pct_OtherNeu")
# Seurat default p_val_adj = Bonferroni x 全 RNA dim. 补充子集 BH (在测试基因里)
mk$p_val_adj_BH <- p.adjust(mk$p_val, method = "BH")
mk <- mk[order(-mk$avg_log2FC), ]
write.table(mk, file.path(OUT, "module2_FindMarkers_Tac1_vs_Others.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(mk, file.path(OUT, "module2_FindMarkers_Tac1_vs_Others.csv"),
          row.names = FALSE)
n_up   <- sum(mk$avg_log2FC > 0.5 & mk$p_val_adj_BH < 0.05, na.rm = TRUE)
n_down <- sum(mk$avg_log2FC < -0.5 & mk$p_val_adj_BH < 0.05, na.rm = TRUE)
cat(sprintf("    Tac1+ 相对 Other_neuron: 上 %d / 下 %d (log2FC>0.5 & BH<0.05)\n",
            n_up, n_down))
tac1_marker_up <- mk$gene_symbol[mk$avg_log2FC > 0.5 & mk$p_val_adj_BH < 0.05]

# 4. 模块 2.2 GO 富集 ------------------------------------------------------
cat(">>> 模块 2.2 Tac1+ 特异 marker GO 富集\n")
do_goana <- function(gene_set, universe_set, label) {
  m_eg <- AnnotationDbi::select(org.Mm.eg.db,
            keys = unique(c(gene_set, universe_set)),
            keytype = "SYMBOL",
            columns = c("ENTREZID","SYMBOL"))
  m_eg <- m_eg[!is.na(m_eg$ENTREZID) & !duplicated(m_eg$SYMBOL), ]
  de_eg  <- unique(m_eg$ENTREZID[m_eg$SYMBOL %in% gene_set])
  uni_eg <- unique(m_eg$ENTREZID[m_eg$SYMBOL %in% universe_set])
  if (length(de_eg) < 5) {
    cat(sprintf("    [%s] DE ENTREZ < 5, 跳过\n", label))
    return(data.frame())
  }
  go <- limma::goana(de = de_eg, universe = uni_eg, species = "Mm")
  go <- go[order(go$P.DE), ]
  go$FDR <- p.adjust(go$P.DE, "BH")
  go$ratio <- go$DE / go$N
  cat(sprintf("    [%s] DE %d, GO terms: %d, FDR<0.05 = %d\n",
              label, length(de_eg), nrow(go), sum(go$FDR < 0.05, na.rm = TRUE)))
  go
}

universe_neu <- mk$gene_symbol  # 所有测试到的基因 = 富集的 universe
go_tac1 <- do_goana(tac1_marker_up, universe_neu, "Tac1+ marker_up")
write.table(go_tac1, file.path(OUT, "module2_GO_Tac1_marker_up.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE)

# 重点关注的咳嗽 / 神经免疫通路关键词
cough_keywords <- "viral|virus|infl|immune|cytokine|interferon|chemokine|pain|sensory|excit|glutamat|synap|cough|nocicept|neuron|axon"
go_focus <- go_tac1[grepl(cough_keywords, go_tac1$Term, ignore.case = TRUE), ]
go_focus <- head(go_focus[go_focus$Ont == "BP", ], 30)
write.table(go_focus, file.path(OUT, "module2_GO_focus_cough_immune.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE)

# 5. 模块 2.3 GSE161878 上调 DEG 在 Tac1+ marker 中 hypergeometric ---------
cat(">>> 模块 2.3 GSE161878 上调 DEG 在 Tac1+ marker 中 enrichment\n")
all_genes_overlap <- intersect(deg_full$SYMBOL, mk$gene_symbol)
N    <- length(all_genes_overlap)
K    <- sum(flu_up %in% all_genes_overlap)            # 流感 up 在测试 universe
n    <- sum(tac1_marker_up %in% all_genes_overlap)    # Tac1+ marker_up 在 universe
k    <- sum(intersect(flu_up, tac1_marker_up) %in% all_genes_overlap)  # 交
hyp_p <- if (n > 0 && K > 0) phyper(k - 1, K, N - K, n, lower.tail = FALSE) else NA
overlap_genes <- intersect(flu_up, tac1_marker_up)
cat(sprintf("    universe=%d, flu_up=%d, Tac1_up=%d, 交=%d, p_hypergeo=%.3g\n",
            N, K, n, k, hyp_p))
cat(sprintf("    交集基因: %s\n",
            paste(head(overlap_genes, 20), collapse=", ")))

# 6. 模块 2.4 核心靶点 VlnPlot --------------------------------------------
cat(">>> 模块 2.4 核心靶点共定位 VlnPlot\n")
core_show <- c("Cx3cr1","Cx3cl1","Gpr151","Gstm1","Klf2","Arrb1",
               "Cd209a","Btla","Il16","Tac1","Slc17a6","Trpv1","Calca","Scn10a")
core_show <- intersect(core_show, rownames(seu_neu))
Idents(seu_neu) <- factor(seu_neu$nlabel, levels = c("Tac1+Vglut2+","Other_neuron"))
p_vln <- VlnPlot(seu_neu, features = core_show, ncol = 4, pt.size = 0.05,
                 cols = c("Tac1+Vglut2+" = "#D62728", "Other_neuron" = "#7F7F7F")) &
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5))
n_panel <- length(core_show); ncol_p <- 4; nrow_p <- ceiling(n_panel/ncol_p)
ggsave(file.path(OUT, "module2_core_targets_violin.pdf"), p_vln,
       width = 3.2 * ncol_p + 1, height = 2.7 * nrow_p + 1,
       device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(OUT, "module2_core_targets_violin.png"), p_vln,
       width = 3.2 * ncol_p + 1, height = 2.7 * nrow_p + 1,
       dpi = 200, type = "cairo", limitsize = FALSE)

# 7. 模块 3 跨数据集通路一致性 -------------------------------------------
cat(">>> 模块 3.1-3.3 GSE161878 vagal up GO + 跨数据集韦恩\n")
go_vagus <- do_goana(flu_up, deg_full$SYMBOL, "GSE161878 vagal up")
write.table(go_vagus, file.path(OUT, "module3_GO_GSE161878_vagal_up.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE)

# 共享 GO BP terms (FDR<0.10 双方)
top_vagus <- rownames(go_vagus)[go_vagus$FDR < 0.10 & go_vagus$Ont == "BP"]
top_tac1  <- rownames(go_tac1) [go_tac1$FDR  < 0.10 & go_tac1$Ont == "BP"]
shared <- intersect(top_vagus, top_tac1)
only_v <- setdiff(top_vagus, top_tac1)
only_t <- setdiff(top_tac1,  top_vagus)
cat(sprintf("    通路 (BP, FDR<0.10): vagus=%d, Tac1+=%d, 共享=%d\n",
            length(top_vagus), length(top_tac1), length(shared)))

# 韦恩 (简易: 用 ggplot 画 2-set)
venn_df <- data.frame(
  label = c(sprintf("仅 vagus\n(%d)", length(only_v)),
            sprintf("共享\n(%d)",       length(shared)),
            sprintf("仅 Tac1+\n(%d)",   length(only_t))),
  x = c(-1, 0, 1), y = 0,
  size = c(length(only_v), length(shared), length(only_t)),
  stringsAsFactors = FALSE
)
if (requireNamespace("ggforce", quietly = TRUE)) {
  library(ggforce)
  p_venn <- ggplot() +
    geom_circle(data = data.frame(x0 = c(-0.5, 0.5),
                                  y0 = c(0, 0), r = c(1, 1),
                                  set = c("Vagus","Tac1+")),
                aes(x0 = x0, y0 = y0, r = r, fill = set),
                alpha = 0.3, color = "black",
                inherit.aes = FALSE) +
    geom_text(data = venn_df, aes(x, y, label = label),
              size = 5, fontface = "bold") +
    scale_fill_manual(values = c(Vagus = "#1F77B4", `Tac1+` = "#D62728")) +
    coord_fixed() +
    theme_void(base_size = 11) +
    labs(title = "GO BP 通路一致性 (FDR<0.10): GSE161878 vagal up ↔ Tac1+ NTS marker_up",
         subtitle = sprintf("共享 %d 通路", length(shared)))
} else {
  # fallback 文字框
  p_venn <- ggplot(venn_df, aes(x, y)) +
    geom_label(aes(label = label), size = 5, fontface = "bold",
               fill = c("#1F77B4","#9467BD","#D62728"),
               alpha = 0.3, color = "white") +
    xlim(-2, 2) + ylim(-1, 1) +
    labs(title = "GO BP 通路一致性 (简化展示)",
         subtitle = sprintf("vagus only=%d, 共享=%d, Tac1+ only=%d",
                            length(only_v), length(shared), length(only_t))) +
    theme_void(base_size = 11)
}
ggsave(file.path(OUT, "module3_pathway_venn.pdf"), p_venn,
       width = 8, height = 4, device = cairo_pdf)
ggsave(file.path(OUT, "module3_pathway_venn.png"), p_venn,
       width = 8, height = 4, dpi = 200, type = "cairo")

# 共享 GO 列表 + 热图
if (length(shared) > 0) {
  shared_df <- data.frame(
    GO_id = shared,
    Term = go_vagus$Term[match(shared, rownames(go_vagus))],
    Ont  = go_vagus$Ont [match(shared, rownames(go_vagus))],
    FDR_vagus = go_vagus$FDR[match(shared, rownames(go_vagus))],
    FDR_Tac1  = go_tac1$FDR[match(shared, rownames(go_tac1))],
    DE_vagus  = go_vagus$DE[match(shared, rownames(go_vagus))],
    DE_Tac1   = go_tac1$DE[match(shared, rownames(go_tac1))],
    stringsAsFactors = FALSE
  )
  shared_df <- shared_df[order(pmax(shared_df$FDR_vagus, shared_df$FDR_Tac1)), ]
  write.table(shared_df, file.path(OUT, "module3_shared_pathways.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  # 热图: 共享 top 25
  top25 <- head(shared_df, 25)
  hm_long <- data.frame(
    Term = rep(top25$Term, 2),
    dataset = rep(c("Vagus_GSE161878","Tac1+_NTS"), each = nrow(top25)),
    negLog10FDR = c(-log10(top25$FDR_vagus + 1e-10),
                    -log10(top25$FDR_Tac1  + 1e-10))
  )
  hm_long$Term <- factor(hm_long$Term, levels = rev(top25$Term))
  p_hm <- ggplot(hm_long, aes(dataset, Term, fill = negLog10FDR)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.1f", negLog10FDR)), size = 3) +
    scale_fill_gradient(low = "white", high = "#D62728", name = "-log10(FDR)") +
    labs(x = NULL, y = NULL,
         title = "Top 25 共享 GO BP 通路: Vagus ↔ Tac1+ NTS",
         subtitle = "颜色 = -log10(FDR)") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 8))
  ggsave(file.path(OUT, "module3_shared_pathway_heatmap.pdf"), p_hm,
         width = 9, height = max(5, nrow(top25) * 0.3 + 2),
         device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(OUT, "module3_shared_pathway_heatmap.png"), p_hm,
         width = 9, height = max(5, nrow(top25) * 0.3 + 2),
         dpi = 200, type = "cairo", limitsize = FALSE)
}

# 8. 模块 3.4 跨数据集核心信号传导分子清单 ---------------------------------
cat(">>> 模块 3.4 跨组织信号传导分子\n")
shared_signal_genes <- intersect(flu_up, tac1_marker_up)
shared_genes_tbl <- data.frame(
  gene_symbol = shared_signal_genes,
  vagus_log2FC = deg_full$log2FC_shrink[match(shared_signal_genes, deg_full$SYMBOL)],
  vagus_FDR    = deg_full$FDR         [match(shared_signal_genes, deg_full$SYMBOL)],
  Tac1_avg_log2FC = mk$avg_log2FC     [match(shared_signal_genes, mk$gene_symbol)],
  Tac1_p_adj_BH   = mk$p_val_adj_BH   [match(shared_signal_genes, mk$gene_symbol)],
  Tac1_pct        = mk$pct_Tac1       [match(shared_signal_genes, mk$gene_symbol)],
  Other_pct       = mk$pct_OtherNeu   [match(shared_signal_genes, mk$gene_symbol)],
  stringsAsFactors = FALSE
)
shared_genes_tbl <- shared_genes_tbl[order(-shared_genes_tbl$vagus_log2FC), ]
write.table(shared_genes_tbl, file.path(OUT, "module3_shared_signaling_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 9. 摘要 ------------------------------------------------------------------
sink(file.path(OUT, "module2_3_summary.txt"))
cat("=== 模块 2/3: Tac1+ 流感富集 + 跨数据集通路一致性 ===\n")
cat(sprintf("时间: %s\n", TS))
cat("脚本: scripts/18_module2_3_enrichment.R\n\n")

cat("=== 模块 2 ===\n")
cat(sprintf("Tac1+Vglut2+ vs Other_neuron FindMarkers (Wilcoxon):\n"))
cat(sprintf("  上 %d / 下 %d (log2FC>0.5 & BH<0.05)\n", n_up, n_down))
cat("\n  Top 20 Tac1+ marker_up:\n")
print(head(mk[mk$avg_log2FC > 0.5 & mk$p_val_adj_BH < 0.05,
              c("gene_symbol","avg_log2FC","p_val_adj_BH",
                "pct_Tac1","pct_OtherNeu")], 20),
      row.names = FALSE)

cat("\n  GO 富集 (Tac1+ marker_up, BP):\n")
top_bp <- head(go_tac1[go_tac1$Ont == "BP", ], 20)
print(top_bp[, c("Term","DE","N","ratio","P.DE","FDR")], row.names = FALSE)

cat("\n  Focus 通路 (cough/immune/sensory/glu/synap):\n")
print(head(go_focus[, c("Term","DE","N","ratio","P.DE","FDR")], 20),
      row.names = FALSE)

cat(sprintf("\n  GSE161878 上调 DEG 在 Tac1+ marker 中 hypergeometric:\n"))
cat(sprintf("    universe=%d, flu_up=%d, Tac1_up=%d, 交=%d, p=%.3g\n",
            N, K, n, k, hyp_p))
cat(sprintf("    重叠基因: %s\n", paste(overlap_genes, collapse=", ")))

cat("\n=== 模块 3 ===\n")
cat(sprintf("GSE161878 vagal up GO BP 富集 (FDR<0.10): %d\n", length(top_vagus)))
cat(sprintf("Tac1+ marker_up GO BP 富集 (FDR<0.10): %d\n", length(top_tac1)))
cat(sprintf("共享通路: %d (vagus only %d, Tac1+ only %d)\n",
            length(shared), length(only_v), length(only_t)))

cat("\nTop 20 共享通路:\n")
if (exists("shared_df") && nrow(shared_df) > 0) {
  print(head(shared_df, 20), row.names = FALSE)
}

cat(sprintf("\n跨组织信号分子清单 (vagus up ∩ Tac1+ marker_up): %d 基因\n",
            nrow(shared_genes_tbl)))
print(shared_genes_tbl, row.names = FALSE)
sink()
cat("\n", paste(readLines(file.path(OUT, "module2_3_summary.txt")), collapse="\n"), sep="")

saveRDS(list(mk = mk, go_tac1 = go_tac1, go_vagus = go_vagus,
             tac1_marker_up = tac1_marker_up,
             shared_pathways = if (exists("shared_df")) shared_df else NULL,
             shared_signal_genes = shared_genes_tbl,
             hyp_test = list(N = N, K = K, n = n, k = k, p = hyp_p,
                             overlap = overlap_genes)),
        file.path(OUT, "module2_3_results.rds"))

cat(sprintf("\n>>> 输出: %s\n", OUT))
