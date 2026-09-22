# 丙氨酸结构交付索引

**请先阅读 [REPORT_zh.md](REPORT_zh.md)。最终推荐模型是本目录的 `alanine.cif`，不是任何 `solution*` 或失败重跑目录中的模型。**

## 最终结构与反射

| 文件 | 用途 |
|---|---|
| `alanine.cif` | 最终候选结构：晶胞、对称、坐标、ADP、统计和限制 |
| `alanine.fcf` | 657个允许反射的Fo²、σ、Fc²及相位，观测值已除以精修比例因子 |
| `observed.hkl` | h k l I σ通用文本，含结尾全零记录；无软件调用含义 |
| `alanine_molecule.xyz` | 一个完整非对称单元分子的笛卡尔坐标，Å；不包含晶胞周期信息 |
| `refined_model.json` | 不舍入的完整模型、协方差参数化说明和统计 |
| `refinement_covariance.npy` | 57参数协方差；顺序见 `refine.py` |
| `refinement_reflections.csv` | 全部精修反射、权重及残差；包括弱/负反射 |
| `bond_geometry.json` | 键长与条件su |
| `difference_density.npy` / `difference_statistics.json` | 有限实测反射集上的差值Fourier |

## 可审计的数据还原

- `frame_metadata.json`：616帧角度、计数诊断、SHA-256。
- `peaks.npy`、`indexed_peaks.npy`、`calibration_centroids.csv`、`geometry.json`：观测峰及几何校准。
- `integrated_run*.npy`、`profiles_run*.npy`：未合并积分和逐帧摇摆记录。
- `unmerged_integrated.csv`：所有几何候选及原始净计数；包含后来未采用的观测。
- `static_masks.npy`、`mask_settings.json`、`pixel_audit_run*.npz`：从原始背景建立的遮挡掩膜及依据。
- `observation_flags.csv`：每个原始积分行的排除理由位标记。
- `masked_observations.csv`、`rejected_observations.csv`：遮挡和仅基于对称等价的异常观测，分别保存。
- `scaled_observations.csv`、`merged_all.csv`、`merging_statistics.json`：比例校正、合并、系统消光和完整度。
- 所有NPY列布局见对应的自编Python文件；CSV文件首行给出列名。

## 检验与重现

- `validation.json`：结构、对称、CIF/FCF读回、接触、收敛及五折条件预测。
- `geometry_validation.json`：晶胞su、几何雅可比秩及留出峰检验。
- `algorithm_tests.json`：压缩、孔径、几何逆变换、输入哈希自测。
- `sensitivity_summary.json`、`sensitivity/`：更宽角窗、更大孔径的独立重积分及结果。
- **`reproduction_verified/`**：已成功执行的全新原始帧→结构流程。`logs/`为每一步stdout/stderr；`run_provenance.json`为解释器、包版本、代码哈希和运行时间。
- `reproducibility_comparison.json`：主模型和独立重求解的一致性检验，允许合法的对称/原点/手性等价变化。
- `MANIFEST.sha256`、`delivery_manifest.json`：本交付文件校验清单（不含清单自身与Python缓存）。

在任务根目录运行：

```powershell
& "H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921\toolbox\pyplain\Scripts\python.exe" "deliverables\run_pipeline.py" --output "deliverables\rerun_new"
```

输出子目录必须尚不存在；驱动程序拒绝覆盖既有结果。`inputs/`应保留在任务根目录，脚本保留在 `deliverables/`。重跑不读取主模型坐标作为起点，而是重新进行分子全局搜索。约需1–2分钟，不安装、不联网。

## 自编主程序

1. `raw_frames.py`：TY6解码与峰提取。
2. `geometry_probe.py`、`geometry.py`：几何约定确认与晶格/模块校准。
3. `audit_pixels.py`、`integrate.py`：原始背景掩膜、Ewald预测、计数积分。
4. `merge.py`：期望计数方差、等价比例校正、合并和消光。
5. `search_molecule.py`：理想连接模板的多起点直接空间求解。
6. `refine.py`：C/N指认比较、F²全矩阵精修、骑乘氢、ADP和差值图。
7. `export_structure.py`、`validate_geometry.py`、`validate_model.py`、`test_algorithms.py`：导出与自检。
8. `run_pipeline.py`、`run_sensitivity.py`、`compare_reproduction.py`：复现与敏感性驱动。

`solve.py`中的显式散射因子和对称常数被主流程复用；该文件的旧电荷翻转入口不属于已验证主流程。所有程序均为本轮自行编写；没有使用晶体学程序库。

## 保留但未采用的探索结果

- `solution.json`、`solution_density.npy`、`solution_trials.json`：早期**未收敛的电荷翻转**，不是最终结构。
- `reproduction/`：首次完整重现暴露的**C/N指认局部极小值**，R1约8%，不是最终结构。该问题已在当前代码和 `reproduction_verified/` 中解决。
- `baseline_before_pixel_mask.*`、`history/`：背景遮挡与计数权重改进前的检查点。
- `decode_probe.py`、`decode_rows.py`、`probe_inputs.py`和各 `audit_*.py`：格式探索/诊断代码；其中中间假设不等同于最终算法。
- `molecular_solution.json`：求解初值，不是最终精修坐标；主目录初值来自早期较小的强反射子集，而完整重现会重新搜索最终数据。

绝对构型未知；尚有显著强度残差、高分辨壳不完整、缺失实验信息及未做独立标准验证等限制。详见报告，不应只凭R因子认定已经可以直接投稿。
