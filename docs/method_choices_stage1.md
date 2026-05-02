# Stage 1 方法选型与可行性评估

> 数据集:GSE42639 (Illumina Mouse v2 芯片,quantile normalized,基因中位数 >50 过滤后保留 N 探针)

## 1. 轻症 / 重症的定义

### 候选与决策

| 候选标签 | 来源 | DEG 数(Mo vs Sham) | 备注 |
|---|---|---|---|
| TX91 | 低毒力 A/Texas/1991 株 | **147**(↑111 / ↓36) | 信号弱但真实,生物学上确为轻症 |
| 0.6LD50PR8 | 0.6 倍 LD50 PR8 株 | 2,012 | "低剂量"但仍是高致病株,信号已很强 |
| 10LD50PR8 | 10 倍 LD50 | 2,981 | 重症 |
| 100LD50PR8 | 100 倍 LD50 | **3,394** | 最重症,信号最强 |

### 决策
- **轻症 = TX91**
  - 理由:(a) 报告 02 中已标注为"轻症流感";(b) DEG 数仅 147,符合"宿主免疫应答有限"的轻症生物学;(c) 不同毒株(TX91 vs PR8)避免了"剂量效应"和"保护机制"混淆。
- **重症 = 100LD50PR8**
  - 理由:(a) Mo 细胞 DEG 最多(3,394),信号最强;(b) 与 TX91 在毒株、剂量两个维度都拉开最大对比度,最容易筛出"轻症特异 / 重症缺失"基因。

### 备选方案(若结果不理想)
- 重症改用 10LD50PR8(剂量更接近,可看是否仍能区分)
- 轻症改用 0.6LD50PR8(同一毒株不同剂量;但与重症共享 PR8 应答主轴,可能稀释保护信号)

---

## 2. 差异分析方法

### 选 limma 的理由
- 数据是芯片(Illumina v2),log 正态近似已成立
- 每组 n=3,小样本:limma 的 empirical Bayes 方差收缩比 t-test、edgeR exact 都更稳
- 已在 `step1-3_GSE42639_analysis_by_celltype.R` 中验证可用

### 设计:联合模型 vs 独立两次拟合
**采用联合模型**(单次 lmFit 在 Mo 全部 15 样本上):
```r
group <- factor(targets_Mo$Treatment,
                levels = c("Sham","TX91","0.6LD50PR8","10LD50PR8","100LD50PR8"))
design <- model.matrix(~ 0 + group)
fit <- lmFit(mat_Mo, design)
contrasts <- makeContrasts(
  C1_TX91_vs_Sham      = groupTX91       - groupSham,
  C2_TX91_vs_100LD50   = groupTX91       - group100LD50PR8,
  levels = design
)
fit2 <- contrasts.fit(fit, contrasts) |> eBayes()
```

**为什么联合不分离**
- 方差估计跨 15 样本汇总,自由度更高(独立两次只有 4 自由度,联合模型 10)
- 两个 contrast 的 p 值在同一参考分布下可比

---

## 3. 显著性阈值

### n=3 的现实
- adj.P.Val < 0.05 + |logFC|≥1 在 n=3 严苛设计下可能漏掉边缘真信号
- 但 TX91 vs Sham 已在严苛阈值下产出 147 个,说明阈值可用

### 主阈值
- adj.P.Val < 0.05 且 logFC ≥ 1(上调)

### 敏感性分析(同时输出)
- 宽松:P.Value < 0.05 且 logFC ≥ 0.585(2^0.585 ≈ 1.5 倍)
- 严苛:adj.P.Val < 0.05 且 logFC ≥ 2

最终交集 A 取主阈值;在报告中附敏感性分析的交集大小作为参考。

---

## 4. 交集 A 的取法

```
A = { gene |
        gene 在 C1 (TX91 vs Sham) 中 adj.P<0.05 且 logFC ≥ 1
      AND gene 在 C2 (TX91 vs 100LD50PR8) 中 adj.P<0.05 且 logFC ≥ 1
}
```

**要求两个 contrast 都 logFC > 0** — 即 TX91 同时高于 Sham 和高于重症。

> ❗ 不能仅看 |logFC| — 必须同向。否则会混入"重症超高、轻症中等、Sham 低"这种非保护型基因。

---

## 5. 风险与回退

| 风险 | 回退方案 |
|---|---|
| 交集 A < 20 基因(太少不足以做下游) | 放宽到敏感性阈值;或重症改用 10LD50PR8 |
| 交集 A > 500(过滤无力) | 加严 logFC ≥ 2;或要求 TX91 logFC 高于 100LD50PR8 logFC ≥ 1 |
| 探针-基因映射缺失严重 | 用 `illuminaMousev2.db` 的最新版本,记录映射前后基因数 |
| Mo 在 100LD50PR8 受损过重导致表达整体偏移 | 在 PCA 上检查 Mo 5 组是否大致连续,若 100LD50 离群严重则换 10LD50PR8 |

---

## 6. 自检清单(脚本运行前)

- [ ] Mo 子集样本数确认 = 15(Sham 3, TX91 3, 0.6 3, 10 3, 100 3)
- [ ] PCA 检查 Mo 5 组分离合理
- [ ] 探针注释包加载成功,映射率 ≥ 60%
- [ ] 输出目录 `results/GSE42639_Mo_protective/` 存在
- [ ] set.seed 已固定(用于任何随机化步骤)
