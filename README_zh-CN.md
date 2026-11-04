# 高阶共现分析与 SMETANA

从多个样本的丰度数据识别显著共现的物种组合，再使用代谢模型评估这些组合的代谢互作潜力。

[English README](README.md) · [安装](INSTALL.md) · [方法与局限](docs/METHODS.md)

运行流程：**样本质控与分组 → 丰度阈值定义存在/缺失 → 跨样本高阶共现推断 → 匹配代谢模型 → SMETANA**。`Computing co-occurrence_Single.py` 根据 [Machado 等（2021）论文](https://doi.org/10.1038/s41559-020-01353-4)的迭代搜索策略处理一份多样本丰度表；`Smetana_FWQ.py` 使用一致的 SBML 模型库模拟入选群落。共现表示统计关联，代谢模型用于估计互作潜力。

## 输入格式与单位

第一行为表头，第一列为唯一物种 ID，其余列为不同样本，物种按行排列：

```text
species_id	sample_01	sample_02	sample_03
species_A	100	0	50
species_B	120	0	60
species_C	0	200	90
```

三样本表仅说明格式，不能满足默认至少 10 个共现样本的要求。数据必须完整、非负、有限且为数值，物种名和样本名不得重复或为空。默认 TSV；CSV 加 `--sep comma`。

- `--input-kind counts`：默认，每个样本按**完整输入列总数**转为相对丰度，再应用阈值。应提供分母完整的计数表。
- `--input-kind tpm`：TPM 除以 1,000,000，不会将筛选过的物种子表重新归一化。列总和不能超过 1,000,000，仅允许小幅舍入误差。
- `--input-kind relative`：输入 0–1 相对丰度，每列总和不超过 1，不会重新归一化子表。百分数需先除以 100。
- `--input-kind presence-absence`：输入严格的 0/1 存在缺失表，丰度阈值由上游定义。

相对丰度**严格大于**阈值才算存在。TPM 数据应明确使用 `--input-kind tpm`，保持物种子表的原始丰度分母。所有样本列均保留，包括全零样本；失败样本需先按明确的质控规则剔除，不能直接当成生物学缺失。

## 安装环境

使用 Python 3.10 或更高版本。通过提供的 Conda 环境安装统计推断和 SMETANA 所需依赖：

```bash
conda env create -f environment_smetana.yml
conda activate smetana_env
```

代谢模拟需要配置受支持的求解器，见 [INSTALL.md](INSTALL.md)。如果只做统计推断，可在已有 Python 环境中运行 `python -m pip install -r requirements.txt`。

## 用 48 个样本的人工 TPM 数据验证

安装见 [INSTALL.md](INSTALL.md)。在仓库根目录运行：

```bash
python "Computing co-occurrence_Single.py" "examples/synthetic_tpm_48.tsv" --input-kind tpm --output "results/demo"
```

该人工数据预期得到 A/B、A/C、B/C 和 A/B/C 四个组合；它们均在 48 个样本中的 12 个样本共现。物种对的期望共现次数为 3，三物种组合为 0.75。该示例只验证计算行为，不代表真实生物学结论；仓库没有可供真实模拟的 SBML 模型。

实际数据建议组织为：

```text
data/
  abundance_tpm.tsv
  models/
    species_A.xml
    species_B.xml
    ...
```

```bash
python "Computing co-occurrence_Single.py" "data/abundance_tpm.tsv" --input-kind tpm --models-dir "data/models" --abundance-cutoff 0.001 --max-combo-size 5 --top-n 1000 --population-size 10000 --seed 42 --output "results/tpm_0p1"
```

`--models-dir` 可省略。提供时，只分析物种 ID 与模型文件名去掉小写 `.xml` 后完全一致的物种；丰度单位转换和阈值处理在模型筛选前完成。该筛选会改变检验的物种集合和多重检验范围，排除情况记录在摘要中。省略时检验全部输入物种，模拟前需确保输出的每个物种均有合适模型。

每次使用新的或空的输出目录；省略 `--output` 时默认目录为 `cooccurrence_results`。`0.001` 是 0.1%，对应**大于 1,000 TPM**；可在另一个目录用 `0.0001`（0.01%，大于 100 TPM）作敏感性分析。这两个阈值不会自动将群落判定为竞争型或合作型。

## 统计与搜索步骤

1. 枚举不重复物种对，统计所有成员共同存在的样本数。
2. 在独立出现的零假设下，将各物种出现概率相乘，乘以样本数得到期望共现次数。
3. 计算单侧二项富集检验，对**本阶全部已评估组合**作 Benjamini-Hochberg 校正，再进行支持度与富集倍数筛选。
4. 同时满足“至少 10 个样本共现、观察值至少为期望值的 2 倍、校正后 p 值不超过 0.05”的组合作为显著结果。48 个样本时，默认支持度意味着至少在 20.8% 的样本中出现。
5. 每阶按共现次数保留最多 1,000 个结果；从本阶显著组合中，按共现次数加权、不放回抽取最多 10,000 个用于扩展。每个入选组合与全部候选物种逐一扩展，去重后评估下一阶。

可用 `--min-samples`、`--min-fold`、`--fdr-cutoff`、`--top-n`、`--population-size` 设置筛选条件和搜索规模。1,000 是**每阶上限**。默认最大规模为 5，可提高；高阶搜索是启发式搜索，不保证找到全部显著组合。`--max-candidates` 默认 1,000,000，超出时明确报错。候选与统计结果会跨阶保留在内存中，因此这个上限是资源保护，实际内存需求取决于各阶候选数量。

默认二项富集尾概率为 `P(X >= 观察次数)`，扩展抽样权重为共现次数。`--pvalue-tail strict` 使用 `P(X > 观察次数)`，`--selection-weight sqrt-frequency` 使用共现次数平方根权重。完整定义、假设及与 [HiOrCo](https://github.com/cdanielmachado/HiOrCo) 的关系见 [METHODS.md](docs/METHODS.md)。

## 输出与 SMETANA 衔接

- `communities.tsv`：已完成各阶的入选群落，**无表头、两列、制表符分隔**，作为 SMETANA 群落定义。
- `size_N_communities.tsv`：仅第 N 阶入选群落。
- `size_N_statistics.tsv`：该阶全部候选的物种列表、观察/期望次数、富集倍数、p 值、校正后 p 值及筛选标记。
- `run_summary.json`：输入哈希、参数、预处理信息、各阶数量和停止原因。

无显著结果时，输出空群落文件并记录原因，保留指定的筛选标准。SMETANA 运行脚本拒绝空群落文件。统计结果以 TSV/JSON 保存，可用支持 TSV 的工具查看。

针对上述目录结构，用共享模型库模拟推断群落：

```bash
conda activate smetana_env
python Smetana_FWQ.py "data/models" --single --communities "../../results/tpm_0p1/communities.tsv" --solver gurobi --mode global --dry-run
python Smetana_FWQ.py "data/models" --single --communities "../../results/tpm_0p1/communities.tsv" --solver gurobi --mode global
```

`--communities` 的相对路径从 `data/models` 解析，也可使用文件绝对路径。`--single` 指一个模型库，不是只分析一个生态样本。`global` 计算 MIP/MRO；`detailed` 计算供体—受体—代谢物层面的得分。可加 `--conda-env smetana_env` 指定环境。

结果与日志位于模型库目录：`global.tsv` 或 `detailed.tsv`、`smetana_run.log`。重复运行可能覆盖旧输出，要比较的模拟请使用独立模型目录副本。预检查只核对 ID 和输入结构，不验证 SBML 内容、模型生长或许可证。

同一套物种模型和环境设置下，每个入选群落模拟一次即可。需要样本特异模型或培养基时，应明确设计并分别保存输出；运行脚本支持样本文件夹和分组文件夹批处理，参数见 `--help`。

## 验证、局限与引用

```bash
python -m unittest discover -s tests -v
```

测试核对统计公式、搜索、输入、可复现性和命令处理；人工数据与模拟进程不等于生物学验证。重复测量、环境混杂、批次差异和模型覆盖偏差应在研究设计中处理。每阶 BH 校正不构成整个自适应搜索的全局 FDR 保证。若研究“竞争—合作极化”，还需一致的 MIP/MRO 比较及同规模随机对照。

`legacy/` 提供可选的数据预处理脚本，以及用于生成人工对照候选群落的 `random_subcommunities.py`。设计对照时需逐一匹配群落规模与数量。环境文件使用 Python 3.10 / SMETANA 1.2.1，安装成功后应导出实际依赖版本以便复现。代码采用 [MIT License](LICENSE)，外部软件使用各自许可；仓库不包含私人数据或求解器密钥。论文、HiOrCo、CarveMe 和 SMETANA 引用见 [英文说明](README.md#validation-limitations-and-citation)。
