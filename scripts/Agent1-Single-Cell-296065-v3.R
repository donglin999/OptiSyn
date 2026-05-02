# ============ 修复版：GSE296065 细胞特异性验证程序 ============
# 修复了Seurat v5的API兼容性问题
# ==========================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(stringr)
library(ggrepel)

# 参数设置
RDS_FILE_PATH <- "C:/docs/2026-04/science/GSE296065/GSE296065_almanzar_flu.rds"
TARGETS_CSV_PATH <- "C:/docs/2026-04/science/GSE296065/core_targets_strict.csv"
OUTPUT_DIR <- "C:/docs/2026-04/science/Specificity_Results_Mo_Macs_v5"

# 筛选标准
FC_THRESHOLD <- 0.5
P_ADJ_THRESHOLD <- 0.05

# 创建输出目录
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("\n", strrep("=", 60), "\n")
cat("     GSE296065 Mo/Mo-Macs特异性验证（Seurat v5兼容版）\n")
cat(strrep("=", 60), "\n\n")

# 1. 加载数据
cat(">>> 步骤1: 加载数据 <<<\n")
if (!file.exists(RDS_FILE_PATH)) {
  stop("❌ 找不到RDS文件:", RDS_FILE_PATH)
}
seurat_obj <- readRDS(RDS_FILE_PATH)
cat("  ✓ 数据加载成功\n")
cat("  细胞数:", ncol(seurat_obj), "\n")
cat("  基因数:", nrow(seurat_obj), "\n")

# 检查Seurat版本
#cat("  Seurat版本:", packageVersion("Seurat"), "\n")
#cat("  SeuratObject版本:", packageVersion("SeuratObject"), "\n")

# 检查默认assay
cat("  默认assay:", DefaultAssay(seurat_obj), "\n")
cat("  可用assay:", Assays(seurat_obj), "\n")

# 检查可用的layer
if ("RNA" %in% Assays(seurat_obj)) {
  cat("  RNA assay的layers:", Layers(seurat_obj[["RNA"]]), "\n")
}

# 2. 识别细胞类型列
cat("\n>>> 步骤2: 识别细胞类型信息 <<<\n")

# 检查元数据结构
cat("元数据列名:\n")
meta_cols <- colnames(seurat_obj@meta.data)
for (i in 1:min(20, length(meta_cols))) {
  cat(sprintf("  %2d. %s\n", i, meta_cols[i]))
}
if (length(meta_cols) > 20) {
  cat(sprintf("  ... 以及另外%d列\n", length(meta_cols) - 20))
}

# 寻找可能的细胞类型列
possible_celltype_cols <- c("celltype", "CellType", "cluster", "seurat_clusters", "ident", 
                           "cell_type", "cell_type_manual", "annotation")

celltype_col <- NULL
for (col in possible_celltype_cols) {
  if (col %in% meta_cols) {
    celltype_col <- col
    cat("  ✓ 找到细胞类型列:", col, "\n")
    break
  }
}

if (is.null(celltype_col)) {
  cat("⚠️ 未找到标准细胞类型列，使用前5列的摘要:\n")
  for (i in 1:min(5, length(meta_cols))) {
    col <- meta_cols[i]
    unique_vals <- unique(seurat_obj@meta.data[[col]])
    cat(sprintf("  %d. %s: %d个唯一值", i, col, length(unique_vals)))
    if (length(unique_vals) <= 5) {
      cat(" (", paste(unique_vals, collapse = ", "), ")", sep = "")
    }
    cat("\n")
  }
  stop("❌ 请手动指定细胞类型列名")
}

# 获取细胞类型信息
seurat_obj$CellType <- as.character(seurat_obj@meta.data[[celltype_col]])

# 显示所有细胞类型
cat("\n所有细胞类型（按数量排序）:\n")
cell_type_summary <- as.data.frame(table(seurat_obj$CellType)) %>%
  arrange(desc(Freq))
colnames(cell_type_summary) <- c("CellType", "Count")

for (i in 1:nrow(cell_type_summary)) {
  cat(sprintf("  %2d. %-30s: %5d 个细胞 (%.1f%%)\n", 
              i, 
              cell_type_summary$CellType[i],
              cell_type_summary$Count[i],
              cell_type_summary$Count[i]/ncol(seurat_obj)*100))
}

# 3. 定义目标群
cat("\n>>> 步骤3: 定义Mo/Mo-Macs相关细胞群 <<<\n")

# 基于您数据中的细胞类型名称，定义Mo/Mo-Macs相关细胞
target_keywords <- c("Monocyte", "Mo", "Macrophage", "Ly6c", "Monocyte.*Ly6c")

# 寻找包含这些关键词的细胞类型
potential_target_groups <- cell_type_summary$CellType[
  grepl(paste(target_keywords, collapse = "|"), cell_type_summary$CellType, ignore.case = TRUE)
]

cat("潜在的目标细胞群（包含单核/巨噬细胞关键词）:\n")
if (length(potential_target_groups) > 0) {
  for (i in 1:length(potential_target_groups)) {
    count <- cell_type_summary$Count[cell_type_summary$CellType == potential_target_groups[i]]
    cat(sprintf("  %2d. %-30s: %5d 个细胞\n", i, potential_target_groups[i], count))
  }
  
  # 自动选择包含"Monocyte"的细胞
  selected_indices <- which(grepl("Monocyte", potential_target_groups, ignore.case = TRUE))
  if (length(selected_indices) == 0) {
    # 如果没有Monocyte，选择第一个包含相关关键词的
    selected_indices <- 1
  }
} else {
  cat("  ⚠️ 未找到包含单核/巨噬细胞关键词的细胞类型\n")
  cat("  请从所有细胞类型中选择:\n")
  for (i in 1:nrow(cell_type_summary)) {
    cat(sprintf("  %2d. %s\n", i, cell_type_summary$CellType[i]))
  }
  selected_input <- readline(prompt = "请输入编号（用逗号分隔）: ")
  selected_indices <- as.numeric(strsplit(selected_input, ",")[[1]])
  potential_target_groups <- cell_type_summary$CellType[selected_indices]
}

# 获取选中的目标细胞类型
target_cell_types <- potential_target_groups[selected_indices]
cat("\n选定的目标细胞群:\n")
for (i in 1:length(target_cell_types)) {
  count <- cell_type_summary$Count[cell_type_summary$CellType == target_cell_types[i]]
  cat(sprintf("  %2d. %-30s: %5d 个细胞\n", i, target_cell_types[i], count))
}

# 背景细胞群：所有其他细胞
background_cell_types <- setdiff(cell_type_summary$CellType, target_cell_types)
cat("\n背景细胞群（其他所有细胞）:\n")
for (i in 1:min(10, length(background_cell_types))) {
  count <- cell_type_summary$Count[cell_type_summary$CellType == background_cell_types[i]]
  cat(sprintf("  %2d. %-30s: %5d 个细胞\n", i, background_cell_types[i], count))
}
if (length(background_cell_types) > 10) {
  cat(sprintf("  ... 以及另外%d种细胞类型\n", length(background_cell_types) - 10))
}

# 4. 读取靶点基因
cat("\n>>> 步骤4: 读取58个靶点 <<<\n")
if (!file.exists(TARGETS_CSV_PATH)) {
  stop("❌ 找不到靶点文件:", TARGETS_CSV_PATH)
}
targets_df <- read.csv(TARGETS_CSV_PATH, stringsAsFactors = FALSE)
target_genes <- unique(targets_df$gene_symbol)
cat("  靶点总数:", length(target_genes), "\n")

# 5. 获取表达数据（Seurat v5兼容方式）
cat("\n>>> 步骤5: 获取表达数据 <<<\n")

# 检查可用的数据层
default_assay <- DefaultAssay(seurat_obj)
cat("  使用assay:", default_assay, "\n")

# 尝试不同的方法获取表达数据
try_get_expression <- function(obj, genes) {
  # 方法1: 使用GetAssayData的layer参数
  if ("data" %in% Layers(obj[[default_assay]])) {
    cat("  使用layer='data'\n")
    return(GetAssayData(obj, layer = "data")[genes, ])
  }
  
  # 方法2: 使用counts层
  if ("counts" %in% Layers(obj[[default_assay]])) {
    cat("  使用layer='counts'\n")
    counts_data <- GetAssayData(obj, layer = "counts")[genes, ]
    # 对数化转换
    return(log1p(counts_data))
  }
  
  # 方法3: 使用旧版slot参数
  tryCatch({
    cat("  尝试旧版API\n")
    return(GetAssayData(obj, slot = "data")[genes, ])
  }, error = function(e) {
    # 方法4: 直接访问assay
    cat("  直接访问assay\n")
    return(obj[[default_assay]]$data[genes, ])
  })
}

# 检查基因是否存在
all_genes <- rownames(seurat_obj)
valid_genes <- intersect(target_genes, all_genes)
missing_genes <- setdiff(target_genes, all_genes)

cat("  存在于数据中的靶点:", length(valid_genes), "/", length(target_genes), "\n")
if (length(missing_genes) > 0) {
  cat("  缺失的靶点（前10个）:\n")
  for (i in 1:min(10, length(missing_genes))) {
    cat(sprintf("    %2d. %s\n", i, missing_genes[i]))
  }
  if (length(missing_genes) > 10) {
    cat(sprintf("    ... 以及另外%d个基因\n", length(missing_genes) - 10))
  }
}

if (length(valid_genes) < 10) {
  cat("⚠️ 警告：可用的靶点基因太少，可能影响分析结果\n")
}

# 获取表达矩阵
cat("  提取表达矩阵...\n")
exp_matrix <- try_get_expression(seurat_obj, valid_genes)
cat("  表达矩阵维度:", dim(exp_matrix), "\n")

# 6. 执行Wilcoxon秩和检验
cat("\n>>> 步骤6: 执行Wilcoxon秩和检验 <<<\n")

# 获取目标细胞和背景细胞的barcodes
target_cells <- colnames(seurat_obj)[seurat_obj$CellType %in% target_cell_types]
background_cells <- colnames(seurat_obj)[seurat_obj$CellType %in% background_cell_types]

cat("  目标群细胞数:", length(target_cells), "\n")
cat("  背景群细胞数:", length(background_cells), "\n")

if (length(target_cells) < 3 || length(background_cells) < 3) {
  stop("❌ 某一组细胞数过少，无法进行有效的统计检验")
}

# 为每个基因计算统计量
results_list <- list()
cat("  正在计算统计量...\n")
pb <- txtProgressBar(min = 0, max = length(valid_genes), style = 3)

for (i in 1:length(valid_genes)) {
  gene <- valid_genes[i]
  
  # 提取表达量
  target_expr <- as.numeric(exp_matrix[gene, target_cells])
  background_expr <- as.numeric(exp_matrix[gene, background_cells])
  
  # 跳过表达量全为0的基因
  if (all(target_expr == 0) && all(background_expr == 0)) {
    setTxtProgressBar(pb, i)
    next
  }
  
  # 执行Wilcoxon检验
  test_result <- tryCatch({
    wilcox.test(target_expr, background_expr, alternative = "greater")
  }, error = function(e) {
    return(list(p.value = NA, statistic = NA))
  })
  
  # 计算平均表达和log2FC
  mean_target <- mean(target_expr, na.rm = TRUE)
  mean_bg <- mean(background_expr, na.rm = TRUE)
  
  # 避免除0错误
  if (mean_bg == 0 || is.na(mean_bg)) {
    log2fc <- ifelse(mean_target > 0, Inf, NA)
  } else {
    log2fc <- log2((mean_target + 1e-10) / (mean_bg + 1e-10))
  }
  
  # 计算表达百分比
  pct_target <- sum(target_expr > 0, na.rm = TRUE) / length(target_expr) * 100
  pct_bg <- sum(background_expr > 0, na.rm = TRUE) / length(background_expr) * 100
  
  results_list[[gene]] <- data.frame(
    gene_symbol = gene,
    mean_target = mean_target,
    mean_background = mean_bg,
    pct_target = pct_target,
    pct_background = pct_bg,
    log2FC = log2fc,
    p_value = test_result$p.value,
    statistic = test_result$statistic,
    stringsAsFactors = FALSE
  )
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 合并结果
if (length(results_list) == 0) {
  stop("❌ 没有有效的基因可进行分析！")
}

stats_df <- bind_rows(results_list)

# 多重检验校正
stats_df <- stats_df %>%
  mutate(
    p_adjust = p.adjust(p_value, method = "BH"),
    is_specific = log2FC > FC_THRESHOLD & p_adjust < P_ADJ_THRESHOLD & !is.infinite(log2FC)
  ) %>%
  arrange(p_adjust)

# 7. 筛选特异性靶点
specific_genes_df <- stats_df %>%
  filter(is_specific) %>%
  arrange(desc(log2FC))

cat("\n>>> 步骤7: 结果统计 <<<\n")
cat("  分析的总基因数:", nrow(stats_df), "\n")
cat("  特异性靶点数 (log2FC >", FC_THRESHOLD, ", p.adjust <", P_ADJ_THRESHOLD, "):", 
    nrow(specific_genes_df), "\n")
cat("  特异性比例:", round(nrow(specific_genes_df)/nrow(stats_df)*100, 1), "%\n")

if (nrow(specific_genes_df) > 0) {
  cat("\nTop 10特异性靶点:\n")
  top_10 <- head(specific_genes_df, 10)
  for (i in 1:nrow(top_10)) {
    cat(sprintf("  %2d. %-15s: log2FC=%.2f, p.adj=%.2e, pct.target=%.1f%%, pct.bg=%.1f%%\n",
                i, top_10$gene_symbol[i], top_10$log2FC[i], top_10$p_adjust[i],
                top_10$pct_target[i], top_10$pct_bg[i]))
  }
} else {
  cat("  ⚠️ 没有发现特异性靶点\n")
}

# 8. 保存结果
write.csv(stats_df, file.path(OUTPUT_DIR, "all_genes_statistics.csv"), row.names = FALSE)
write.csv(specific_genes_df, file.path(OUTPUT_DIR, "specific_targets.csv"), row.names = FALSE)

# 保存靶点存在性信息
existence_stats <- data.frame(
  gene = target_genes,
  exists_in_data = target_genes %in% valid_genes,
  stringsAsFactors = FALSE
)
write.csv(existence_stats, file.path(OUTPUT_DIR, "gene_existence_stats.csv"), row.names = FALSE)

# 9. 可视化
if (nrow(stats_df) > 0) {
  cat("\n>>> 步骤8: 生成可视化图表 <<<\n")
  
  # 9.1 火山图
  volcano_data <- stats_df %>%
    mutate(
      log10_padj = -log10(p_adjust + 1e-300),  # 防止-log10(0)产生Inf
      significance = case_when(
        is_specific ~ "特异性高表达",
        p_adjust < 0.05 & log2FC > 0 ~ "显著但未达阈值",
        p_adjust < 0.05 ~ "显著低表达",
        TRUE ~ "不显著"
      ),
      # 计算折叠变化倍数
      fold_change = round(2^log2FC, 1)
    )

  volcano_plot <- ggplot(volcano_data, aes(x = log2FC, y = -log10(p_adjust), color = significance)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c(
      "特异性高表达" = "#D62728",
      "显著但未达阈值" = "#FF7F0E",
      "显著低表达" = "#1F77B4",
      "不显著" = "gray70"
    )) +
    geom_vline(xintercept = FC_THRESHOLD, linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log10(P_ADJ_THRESHOLD), linetype = "dashed", color = "black", alpha = 0.5) +
    labs(
      title = paste0("Mo/Mo-Macs特异性验证火山图 (n=", nrow(specific_genes_df), "/", nrow(stats_df), ")"),
      subtitle = paste0("目标: ", paste(target_cell_types, collapse = ", ")),
      x = "log2 Fold Change (目标群/背景群)",
      y = "-log10(adjusted p-value)",
      color = "显著性"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "right"
    )
  volcano_plot <- volcano_plot +
    scale_y_continuous(
      limits = c(0, 800),  # 直接设置范围
      expand = expansion(mult = c(0, 0.2))  # 在顶部扩展20%的空间
    )

  
  # 添加Top特异性基因的详细信息标签
  if (nrow(specific_genes_df) > 0) {
    # 选择要显示的Top基因（最多10个）
    top_n_genes <- min(10, nrow(specific_genes_df))
    top_label_genes <- specific_genes_df %>%
      arrange(desc(log2FC)) %>%
      head(top_n_genes)
    
    # 准备标签数据
  

    label_data <- volcano_data %>% 
      filter(gene_symbol %in% top_label_genes$gene_symbol) %>%
      filter(is.finite(log10_padj)) %>%
      left_join(top_label_genes %>% select(gene_symbol, pct_target, `pct_background`), 
                by = "gene_symbol") %>%
      mutate(
        # 创建包含详细信息的标签文本
        label_text = paste0(
          gene_symbol, 
          "\nlog2FC=", sprintf("%.2f", log2FC),
          " (", fold_change, "x)",
          "\np.adj=", formatC(p_adjust, format = "e", digits = 1),
           "\n表达率: ", 
            round(pct_target.y, 1), "% vs ",  # 使用 .y 后缀
            round(pct_background.y, 1), "%"   # 使用 .y 后缀
        )
      )
    


    if (nrow(label_data) > 0) {
      # 方法1：使用ggrepel显示详细标签
      volcano_plot <- volcano_plot +
        ggrepel::geom_text_repel(
          data = label_data,
          aes(label = label_text),
          size = 2.8,  # 稍微调小字号
          color = "black",
          fontface = "bold",
          # 布局控制
          box.padding = 2.8,      # 增加标签周围的空白
          point.padding = 0.5,    # 标签与点的距离
          min.segment.length = 0.2,
          segment.size = 0.4,
          segment.color = "gray40",
          segment.alpha = 0.7,
          # 避让算法控制
          max.time = 2,
          max.iter = 10000,
          force = 2.5,            # 增加排斥力
          force_pull = 0.8,       # 增加吸引力
          direction = "y",
          seed = 123,
          # 边界限制
          xlim = c(NA, NA),
          ylim = c(NA, NA),  # 限制标签不超出Y轴顶部的95
          max.overlaps = Inf
        )
      
      # 在图上添加说明文本
      volcano_plot <- volcano_plot +
        annotate("text", 
                 x = min(volcano_data$log2FC, na.rm = TRUE) + 0.5,
                 y = y_limit * 0.95,
                 label = paste0("Top ", top_n_genes, " 特异性基因"),
                 hjust = 0, vjust = 1, size = 3.5, color = "black", fontface = "bold") +
        annotate("text", 
                 x = min(volcano_data$log2FC, na.rm = TRUE) + 0.5,
                 y = y_limit * 0.90,
                 label = paste0("筛选标准: log2FC > ", FC_THRESHOLD, ", p.adj < ", P_ADJ_THRESHOLD),
                 hjust = 0, vjust = 1, size = 3, color = "gray40")

      volcano_plot <- volcano_plot 
    }
  }


  # 保存火山图
  ggsave(
    file.path(OUTPUT_DIR, "specificity_volcano_plot.png"),
    volcano_plot,
    width = 12,
    height = 10,
    dpi = 300
  )
  cat("  ✓ 火山图已保存: specificity_volcano_plot.png\n")
  
  # 9.2 表达百分比条形图
  if (nrow(specific_genes_df) > 0) {
    top_genes_plot <- specific_genes_df %>%
      arrange(desc(log2FC)) %>%
      head(22) %>%
      mutate(gene_symbol = factor(gene_symbol, levels = rev(gene_symbol)))
    
    pct_plot <- ggplot(top_genes_plot, aes(x = gene_symbol)) +
      geom_bar(aes(y = pct_target, fill = "目标群"), stat = "identity", alpha = 0.7, width = 0.6) +
      geom_bar(aes(y = pct_background, fill = "背景群"), stat = "identity", alpha = 0.7, width = 0.6, 
               position = position_nudge(x = 0.3)) +
      scale_fill_manual(values = c("目标群" = "#D62728", "背景群" = "#1F77B4")) +
      coord_flip() +
      labs(
        title = "Top特异性靶点的表达百分比",
        subtitle = paste0("显示前", nrow(top_genes_plot), "个特异性靶点"),
        x = "基因",
        y = "表达细胞百分比 (%)",
        fill = "细胞群"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.text.y = element_text(size = 9),
        legend.position = "right"
      )
    
    # 单独保存条形图
    ggsave(
      file.path(OUTPUT_DIR, "specificity_expression_percentage.png"),
      pct_plot,
      width = 10,
      height = 8,
      dpi = 300
    )
    cat("  ✓ 表达百分比条形图已保存: specificity_expression_percentage.png\n")
  }
}

# 10. 生成报告
cat("\n>>> 步骤9: 生成分析报告 <<<\n")

report_lines <- c(
  "GSE296065 Mo/Mo-Macs特异性验证分析报告",
  paste("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  strrep("=", 60),
  "",
  "1. 分析概述",
  strrep("-", 40),
  paste("   输入文件:", RDS_FILE_PATH),
  paste("   靶点文件:", TARGETS_CSV_PATH),
  paste("   总细胞数:", ncol(seurat_obj)),
  paste("   总基因数:", nrow(seurat_obj)),
  paste("   分析靶点数:", length(valid_genes), "/", length(target_genes)),
  "",
  "2. 细胞分组",
  strrep("-", 40),
  paste("   目标细胞类型:", paste(target_cell_types, collapse = ", ")),
  paste("   目标细胞数:", length(target_cells)),
  paste("   背景细胞类型数:", length(background_cell_types)),
  paste("   背景细胞数:", length(background_cells)),
  "",
  "3. 分析参数",
  strrep("-", 40),
  paste("   log2FC阈值:", FC_THRESHOLD),
  paste("   校正p值阈值:", P_ADJ_THRESHOLD),
  paste("   统计方法: Wilcoxon秩和检验"),
  paste("   多重检验校正: Benjamini-Hochberg (FDR)"),
  "",
  "4. 分析结果",
  strrep("-", 40),
  paste("   分析基因数:", nrow(stats_df)),
  paste("   特异性靶点数:", nrow(specific_genes_df)),
  paste("   特异性比例:", round(nrow(specific_genes_df)/nrow(stats_df)*100, 1), "%")
)

if (nrow(specific_genes_df) > 0) {
  report_lines <- c(report_lines,
                   "",
                   "5. Top特异性靶点",
                   strrep("-", 40))
  
  for (i in 1:min(10, nrow(specific_genes_df))) {
    gene <- specific_genes_df$gene_symbol[i]
    log2fc <- specific_genes_df$log2FC[i]
    p_adj <- specific_genes_df$p_adjust[i]
    pct_target <- specific_genes_df$pct_target[i]
    pct_bg <- specific_genes_df$pct_background[i]
    
    report_lines <- c(report_lines,
                     paste("   ", i, ". ", gene, 
                           " (log2FC=", round(log2fc, 2), 
                           ", p.adj=", formatC(p_adj, format = "e", digits = 2),
                           ", pct.target=", round(pct_target, 1), "%",
                           ", pct.bg=", round(pct_bg, 1), "%)", sep = ""))
  }
}

# 保存报告
writeLines(report_lines, file.path(OUTPUT_DIR, "analysis_report.txt"))
cat("  ✓ 报告已保存\n")

cat("\n" , strrep("=", 60), "\n")
cat("    ✅ 分析完成！\n")
cat(strrep("=", 60), "\n\n")

cat("📁 输出文件已保存至:", OUTPUT_DIR, "\n")
cat("   ├── all_genes_statistics.csv: 所有基因的统计结果\n")
cat("   ├── specific_targets.csv: 特异性靶点列表\n")
cat("   ├── gene_existence_stats.csv: 基因存在性统计\n")
cat("   ├── analysis_report.txt: 分析报告\n")
if (exists("combined_plot")) {
  cat("   └── specificity_analysis_plots.png: 可视化图表\n")
}