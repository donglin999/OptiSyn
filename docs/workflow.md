# OptiSyn 工作流程(更新版)

> 最后更新:2026-04-27

## 核心机制假设
> **轻症流感 → 肺 Mo 细胞激活保护基因 → 分泌配体 → 迷走神经感觉元受体 → 触发保护性咳嗽反射 → 中药复方调控**

⚠️ **关键约束**:小鼠无典型咳嗽反射,**不能**把 Mo 细胞保护基因直接与 GSE161878(神经数据)取交集 — 不同器官。必须通过"配体—受体"做空间转录组学桥接。

---

## 数据集清单

| 用途 | GSE | 说明 |
|---|---|---|
| 主分析(分细胞) | **GSE42639** | 小鼠流感分细胞类型,Mo 核心数据 |
| 全肺验证 | **GSE31022** | 流感急性期全肺转录组 |
| 旧 | GSE63786 | 当前未启用 |
| 神经-感染 | **GSE161878** | 流感感染 d4 迷走感觉神经元单细胞 |
| 神经图谱 | **GSE124312** | 迷走神经 18 亚型基线单细胞图谱 |

数据集分组与样本量参见 `reports/02_GSE42639_分细胞类型分析报告.md`。

### GSE42639 细胞类型注释(权威,以本表为准)

| 缩写 | 全称 |
|---|---|
| Nh | pulmonary CD45neg(肺 CD45⁻ 非免疫细胞:上皮/内皮/基质) |
| Am | alveolar macrophage(肺泡巨噬细胞) |
| Ne | neutrophil(中性粒细胞) |
| Mo | Ly6Chi mononuclear myeloid cell(Ly6Chi 单核髓系细胞) |
| Ly | lymphocyte(BC/TC/NK 合并) |

---

## Stage 1 — Mo 细胞保护基因(交集 A)

数据源:GSE42639;只用 Mo 子集(15 样本,5 处理 × n=3)。

| 比较 | 设计 | 含义 |
|---|---|---|
| 对比 1 | TX91 vs Sham(Mo) | 轻症流感激活的基因 |
| 对比 2 | TX91 vs 100LD50PR8(Mo) | 轻症特异性高表达 / 重症受损或缺失 |

**轻症 / 重症定义**:TX91 = 轻症(低毒力 Texas/91 株,Mo 中 147 DEG);100LD50PR8 = 重症(Mo 中 3,394 DEG,信号最强)。详见 `docs/method_choices_stage1.md`。

**交集 A 定义**:在两个对比中均显著上调的基因(同向 logFC>0)→ "Mo 细胞轻症流感早期特异性保护基因"。

产物:
- `results/GSE42639_Mo_protective/Mo_TX91_vs_Sham.txt`
- `results/GSE42639_Mo_protective/Mo_TX91_vs_100LD50PR8.txt`
- `results/GSE42639_Mo_protective/Mo_protective_setA.txt` ⭐ 交集 A
- `results/GSE42639_Mo_protective/Mo_protective_setA.rds`

---

## Stage 2 — 全肺信号验证(GSE31022,可选)

把交集 A 的基因投影到 GSE31022 全肺早期时间点,看是否同向上调,衡量 Mo 信号在组织微环境中的权重。

产物:
- `results/validation/setA_in_GSE31022.txt`

---

## Stage 3 — 分泌组过滤(免疫 → 神经的桥)

把交集 A ∩ 分泌蛋白库 → **保护性分泌配体集 L**:
- 数据库来源:CellChatDB(Ligands)、SecretomeP、Human Protein Atlas Secretome
- 保留类型:细胞因子、趋化因子、神经肽、脂质代谢酶

产物:
- `results/secretome/protective_ligands_L.txt`
- `results/secretome/database_versions.txt`(可追溯性)

---

## Stage 4 — 迷走神经受体与通讯网络

1. **受体集 R**:GSE161878 中流感早期(d4)显著上调的受体基因
2. **L–R 配对**:CellChat 或 NicheNet,基于 PPI 计算通讯概率(L 来自 Stage 3,R 来自上一步)
3. **神经元亚型定位**:用 GSE124312 单细胞注释,把配对成功的 R 投到 18 亚型上,看是否富集在 TRPV1⁺ C 纤维 / Aδ 机械感受纤维(咳嗽相关)
4. **下游通道验证**:被激活受体在 GSE161878 中是否引起 TRP / P2X 等神经冲动通道上调

产物:
- `results/vagal/upregulated_receptors_R.txt`
- `results/communication/LR_pairs_cellchat.rds`
- `results/communication/LR_pairs_summary.txt`
- `results/core_axis/Mo_to_vagal_axis.txt` ⭐ 核心通讯轴

---

## Stage 5 — OptiSyn 多维 AI 中药筛选

### 5.1 网络拓扑初筛
- 输入:Stage 3 配体 L + Stage 4 受体 R(核心轴)→ TCMSP / HERB 反向映射
- 指标:Overlap rate、Hscore(随机游走)、Z-score(PPI 邻近度,要求 Z<0)

### 5.2 成分精炼
- AlphaFold 结构 + AutoDock Vina 高通量对接,保留 ΔG < -5 kcal/mol
- Tanimoto + K-means 聚类去冗余(Silhouette ≥ 0.45,Dunn 指数评估)

### 5.3 GCN 协同预测 + 君臣佐使解析
- 节点特征:拓扑(Z) + 化学(对接 + 聚类) + 临床(Mscore / Pscore,来自 CNKI / 专利)
- 多层 GCN → Sigmoid 输出药对协同概率
- SVD 降维 → 前 25% 为君臣(主要靶向 Mo–神经轴),其余佐使

产物:
- `results/tcm/candidate_herbs_topology.txt`
- `results/tcm/docking_scores.txt`
- `results/tcm/gcn_synergy_pairs.txt`
- `results/tcm/final_formula_recommendation.txt`

---

## Agent 分工

| Stage | 主要 agent | 辅助 agent |
|---|---|---|
| 1 | data-analyst | method-selector, data-preprocessor |
| 2 | data-analyst | result-analyst |
| 3 | result-analyst | third-party-validator(数据库版本) |
| 4 | data-preprocessor(单细胞) → data-analyst | method-selector(CellChat vs NicheNet), result-analyst |
| 5 | method-selector + deep-learning-analyst | result-analyst |

---

## 报告结构对应

- `reports/01_GSE31022_流感急性期分析报告.md`(已有)
- `reports/02_GSE42639_分细胞类型分析报告.md`(已有)
- `reports/04_核心靶点筛选报告.md`(已有)
- `reports/05_整体分析汇总报告.md`(已有)
- `reports/06_Mo神经通讯轴分析报告.md`(待写,Stage 1–4)
- `reports/07_OptiSyn中药复方预测报告.md`(待写,Stage 5)
