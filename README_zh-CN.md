# MetaCoNet：群落组合生成与 SMETANA 批量分析

本仓库提供从物种代谢模型生成候选亚群落、再批量运行 SMETANA 的脚本。[English README](README.md) · [安装说明](INSTALL.md)

## 功能与适用范围

- `Computing co-occurrence_Single.py`：读取一个样本或模块目录中的 `.xml` 模型，按不同成员数量分配名额，随机生成互不重复的候选群落，输出 TSV，可选输出 Excel。
- `Smetana_FWQ.py`：检查模型与群落定义，在各样本文件夹内运行 SMETANA，保存日志并报告失败。
- `legacy/`：保留原始数据预处理辅助脚本，使用前需要修改路径并核对数据格式。

**这里的组合生成脚本不做跨样本共现频率统计、二项检验、FDR 校正或丰度筛选。** 它从已选定的模型集合抽样。候选组合不等同于显著共现关系或已验证的互利关系。生态网络划分、模块选择与丰度筛选需在上游完成。

## 输入要求

每个样本目录存放多个有效的 SBML 文件，例如 `species_A.xml`、`species_B.xml`。模型文件名去掉 `.xml` 后就是物种 ID，必须与群落 TSV 中的 ID 完全一致。使用小写 `.xml` 扩展名；包装脚本不处理压缩模型。

```text
data/
  sample_01/
    species_A.xml
    species_B.xml
    species_C.xml
  sample_02/
    ...
```

仓库不包含原始序列、生物数据、实际模型或求解器许可证。

## 快速使用

按 [INSTALL.md](INSTALL.md) 安装软件并配置求解器后，在仓库根目录运行：

```bash
conda activate smetana_env
python "Computing co-occurrence_Single.py" "data/sample_01" --max-combo-size 5 --total-budget 1000 --seed 42 --output "data/sample_01/communities.tsv" --excel "data/sample_01/communities.xlsx"
```

对每个样本分别生成群落文件。默认抽取 2–5 个成员的组合，总数最多 1000；可用组合不足时输出全部可用组合。名额尽可能在不同规模之间均衡分配；固定模型 ID、参数、随机种子和 Python 版本可复现结果。

SMETANA 使用**无表头、两列、制表符分隔**的 TSV：第一列群落名，第二列物种 ID。Excel 仅用于查看，含表头，不可直接作为 SMETANA 输入。

```text
community1	species_A
community1	species_B
community2	species_A
community2	species_C
```

先检查输入和将要执行的命令：

```bash
python Smetana_FWQ.py "data" --solver gurobi --dry-run
```

检查通过后运行全局分析（MIP/MRO）：

```bash
python Smetana_FWQ.py "data" --solver gurobi --mode global
```

详细分析改为 `--mode detailed`。结果 `global.tsv` 或 `detailed.tsv` 和日志 `smetana_run.log` 位于每个样本目录。重复运行可能覆盖这些文件，需要保留的运行请另存输入目录。输入检查不会验证 SBML 内容或求解器许可证；模型质量与结果解释仍需自行核对。

单个样本直接运行：

```bash
python Smetana_FWQ.py "data/sample_01" --single --solver gurobi
```

原始的“分组/样本”目录结构，如 `data/CK_XML/sample_01/`：

```bash
python Smetana_FWQ.py "data" --groups CK_XML DL_XML SL_XML UL_XML --solver gurobi
```

默认从当前环境调用 `smetana`。也可加 `--conda-env smetana_env` 指定环境；`--communities` 可修改群落文件名。某个样本失败后继续处理其他样本，整个任务最终返回失败状态，详情见日志。完整参数见脚本 `--help`。

## 安装与验证说明

本次发布的包装脚本要求 Python >=3.8；环境文件建议 Python 3.10 / SMETANA 1.2.1。历史本地笔记中的 Python 3.6 / SMETANA 1.2.0 属于旧环境，不是本次包装脚本的运行要求。环境文件尚未进行完整安装与生物模拟验证，也不是依赖锁定文件。使用已有环境或更换求解器后请重新验证。

```bash
python -m unittest discover -s tests -v
```

测试覆盖组合数量、可复现性、TSV 格式、输入检查及批处理命令；外部求解器调用使用模拟测试，不能据此宣称生物学结果得到验证。


