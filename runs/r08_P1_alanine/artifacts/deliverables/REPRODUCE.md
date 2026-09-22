# 复算与文件说明

## 1. 最短的结构复算（已经实际验证）

把本目录的 `alanine.ins`、`alanine.hkl` 复制到一个**新的空目录**，在该目录执行：

```powershell
& "H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921\toolbox\shelx\shelxl.exe" alanine
```

本次使用 SHELXL-2019/3。输出应为 964 条独立反射、68 参数、0 restraints，R1=0.0302、wR2=0.0732、S=1.095。实际隔离复跑记录见 `validation/validation.json` 和完整记录包中的 `work/validation_list4/shelxl_rerun/`。请检查 `.lst/.res` 是否完成，不能只看程序退出码：该 SHELXL 在某些输入错误下也返回 0。

- `.ins` 是产生最终结果的计算输入；`.res` 是收敛输出；`.hkl` 是带符号、未合并的 HKLF4 测量数据。
- 最终 `.fcf` 是 **LIST 4**，保留 964 个 Friedel 独立反射，适合与 CIF 一同验证。
- `provenance/alanine_fourier_LIST6.fcf` 是另外保留的差值图文件，合并 Friedel 后为 596 条记录；不要将其替代投稿用 `.fcf`。
- SHELXL 直接重跑生成的 CIF 不含全部人工核实的实验元数据；交付 CIF 仅据已知事实补全，并保留原生 CIF 和逐字段修改记录。
- 最终结构来源 `work/refine/10_publication`，与 `07_final` 使用相同数据、模型、权重和精修设置；只改 LIST 4、增加 CONF 几何表，并仅列三个明确的 N—H···O 氢键。

## 2. 最终降阶文件

`processing/` 保存最终掩膜、积分、尺度化结果及 MTZ、强度审计和晶胞协方差：

- `masked_input.expt` + `refined_joint_esd.refl`：最终重新积分的输入。
- `geometry_0.mask` 至 `geometry_3.mask`、`mask_definition.json`：四组探测器几何的像素掩膜及明确边界。
- `integrated.expt/refl`、`scaled.expt/refl`：中间/最终 DIALS 结果。
- `scaled_unmerged.mtz`、`scaled_merged.mtz`：两种反射格式，保留 Friedel 信息。
- `intensity_audit.txt/json`：从实际 signed HKL 计算的完整度、冗余度、Rmerge、CC1/2、分辨率壳层和异常差统计。
- `common_cell_covariance.npz/json`：约束参数协方差、展开矩阵 T 和完整协方差。

注意 `.expt` 内含本机绝对路径（原始帧、掩膜等）。在别的目录复算时，应在**副本**中将原运行根目录替换为新工作区根目录；原始帧路径应指向相同、经 SHA-256 核验的 inputs。不要修改原始 inputs。仅用最终 HKL 复算 SHELXL 不依赖这些绝对路径。

## 3. 从原始帧复核处理路线

`calculation_records.zip` 保留原来的 `work/` 层次、成功及失败的控制台/计算文件、所有本地脚本、PHIL、输入哈希和命令清单。解压后各脚本需位于其工作区根目录，因为 `run_command.ROOT` 取脚本所在目录。不要在唯一的存档副本中覆盖重跑。

核心软件为 DIALS 3.30、SHELXT 2018/2、SHELXL 2019/3；解释器为：

```text
C:\users\<user>\miniforge3\envs\dials\python.exe
```

调用 DIALS 时需让其 `Library\bin` 与 `Scripts` 在 PATH 上。`run_command.py` 会为子进程设置该路径并记录命令。工具路径若变化，须在复算副本中调整。以下是**成功路线的顺序说明，而非已测试的一键重放器**；准确参数以 `provenance/commands.jsonl`、PHIL 和源脚本为准。

1. `import_frames.py`：为 646 帧生成逐文件 PHIL 并导入；`prepare_processing.py` 排除 30 个 `pre_` 帧，仅保留六扫描 616 帧并寻峰。
2. `dials.index production.expt strong.refl joint_indexing=True`；`dials.refine_bravais_settings indexed.expt indexed.refl nproc=4`。不需要输入参考晶胞或坐标。选推荐 oP 的 `bravais_setting_5.expt`（本次轴变换为恒等）。
3. 静态面板几何精修：`scan_varying=False detector.panels=hierarchical detector.hierarchy_level=1`；输出 `refined_static.expt/refl`。
4. `refine_and_integrate.py`：为各扫描复制独立取向，固定晶胞做扫描变化精修，完成第一次积分和消光检查。
5. `dials.refine refined.expt refined.refl joint_cell.phil`：精确约束公共正交晶胞；然后在 `work/dials` 运行 `common_cell_covariance.py` 得到同一拟合解及完整 ESD。这一扩展仅在本地进程中生效，不更改安装软件。
6. `process_final_data.py`：未掩膜对照的积分、尺度化、SHELX 导出和 SHELXT。**历史目录名 `work/final_data` 不代表最终接受的数据。**
7. `prepare_refinement.py` 及 `refine_stage.py --add-h`：从独立求解的六个重原子建立初始模型。失败的 `01_nonH`、`02_riding_H` 不使用；成功版本为 `01b_nonH`、`02b_riding_H`。
8. `diagnose_shadows.py`、`evaluate_shadow_background.py`：在原始像素中核查束挡/支撑阴影，按四种几何累积背景。`process_masked_data.py` 用保存的探测器坐标规则生成掩膜，重新积分、尺度化、导出和审计，得到 `work/masked_data`。
9. 用 `refine_stage.py` 将同一模型换用掩膜后的 HKL；`test_hydrogens.py` 做 NH 省略图和独立 N—H 精修；`finalize_refinement.py` 做权重收敛及 0.80 Å/全部实测数据对照。
10. 最终格式用 LIST 4、CONF 和三个明确氢键；`prepare_delivery.py` 记录元数据补全；`validate_structure.py` 记录独立检查。打包/格式脚本含保护现有输出的断言，并有当次处理的先后依赖，不能在存档上盲目顺次覆盖执行。

束挡保护区是根据本次原始像素背景定义的保守掩膜，不是 Fo−Fc 逐条反射剔除。可用保留的未掩膜对照和像素统计复核其依据。所有标准积分/partiality/等价反射离群值过滤均在程序日志中；没有使用 SHELXL 个别反射 OMIT。

## 4. 记录口径与完整性

- `provenance/raw_frames_manifest.json`：646 个原始帧的大小与 SHA-256；任务结束前全部再次核验，未改变。
- `provenance/commands.jsonl`：实际命令、工作目录、起止时间、耗时与退出码。部分科学脚本直接调用的验证子程序另在 `validation.json` 中记录。
- `provenance/cif_metadata_changes.json`：所有 CIF 元数据补全/修正；原坐标、ADP 和内嵌 HKL/RES 不被手工篡改。
- 晶胞长度 ESD 用实际约束精修协方差；体积 ESD 包含一般相关项，为 0.15 Å³。SHELXL 的原生近似值为 0.10 Å³，故本地校验可能提示差异，已有方法说明。
- `SHA256SUMS.json`：最终交付文件校验和（不包含清单本身）。原始帧不重复打包，请与本次 inputs 一同保存。
