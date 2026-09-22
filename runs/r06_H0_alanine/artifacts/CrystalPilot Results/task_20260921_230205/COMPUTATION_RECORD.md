# 计算与复现记录

所有结构生成、编辑、精修和结构交付均经 `mcp__crystalpilot` 工具；没有使用shell直接运行解析/精修引擎。`inputs/alanine`只读。没有读取其他项目、参考结构或评分资料，也没有联网检索该样品答案。

## 原始数据处理

| 顺序 | 工具与参数 | 保存的结果 |
| --- | --- | --- |
| 1 | `import_frames(frames_dir="inputs/alanine/frames", min_run_frames=12, dials_env="C:\users\<user>\miniforge3\envs\dials")` | 6正式扫描，616帧；预实验3x10帧另留输入中 |
| 2 | `find_spots(nproc=1)` | 2821强点，未设置分辨率截断 |
| 3 | `index_frames()` | 无晶胞/空间群先验；P1，2714点，96.2% |
| 4 | `integrate_frames(nproc=1)` | 几何精修及积分；工具n_integrated=4136；表含4422预测记录 |
| 5 | `scale_and_export(composition="C3 H7 N O2")` | 默认combine强度；自动建议P212121；无resolution参数 |
| 6 | `create_start_model(composition="C3 H7 N O2", space_group="P 21 21 21")` | 电荷翻转起模，6非氢位点；节点n0000，数据d000001 |

处理程序DIALS 3.30。缩放日志记录4160观测、687非异常独立反射、Rmerge=0.039、Rmeas=0.042、Rpim=0.015、CC1/2=1.000、平均I/sigma=32.1。这些是处理阶段统计，不代替最终SHELXL统计。最终导出4134行；后续摄入报告2行sigma<=0被自动排除。

原始处理文件完整保留在当前项目 `.crystalpilot/frames/`：`imported.expt`、`strong.refl`、`indexed.expt/.refl`、`refined.expt`、`integrated.expt/.refl`、`scaled.expt/.refl/.mtz`、`dials.ins/.hkl`、`state.json`。工具生成的规范化帧名是输入帧的硬链接，未改像素。`logs/`含import、find_spots、index、refine、integrate、symmetry、scale、export、export_mtz日志。

## 模型链

| 节点 | 操作 |
| --- | --- |
| n0000 | 工具内部电荷翻转粗解；初始C4O2标签不是最终组成 |
| n0001 | `edit_atoms`将原C1重判为N |
| n0002 | `rename_atoms`：原C1->N1、原C3->C1、原C4->C3 |
| n0003 | `set_adp`：6个非氢原子全部各向异性 |
| n0004 | `run_shelxl(mode="adopt", l_s=10)`，尚无H |
| n0005 | `assemble_asu()`，仅合法对称等价拼接 |
| n0006 | `add_hydrogens(elements=["C"], force_kind={"N1":"CH3_rotating","C3":"CH3_rotating"})`，7H |
| n0007 | `run_shelxl(mode="adopt", l_s=15)` |
| n0008 | 输入150 K元数据后，`run_shelxl(mode="adopt_wght", l_s=12, wght_rounds=6)`，实际追加1轮 |
| n0009 | `set_weights(a=0.1393,b=0)`，采纳上轮建议用于核查 |
| n0010 | `run_shelxl(mode="adopt", l_s=12, extra_cards=[3条EQIV,3条HTAB])`，最终匹配作业 |

最终SHELXL输入即 `final.ins`，对应反射为 `final.hkl`，精修输出模型为 `final.res`，结构因子为 `final.fcf`。`final.cif`含SHELXL ACTA几何及统计。需要独立复算时，这一同名INS/HKL组合可直接供SHELXL使用；在本工作台继续修改仍须使用MCP工具，保持节点链。

原始最终作业完整目录：`.crystalpilot/refine/shelxl/job_20260921_231640/`，含`job.ins/.hkl/.res/.lst/.cif/.fcf`及作业记录。SHELXL版本2019/3；最终命令记录为`job -a50000 -b3000 -c624 -g0 -m0 -t4`。这里仅引用已执行的工具日志，不是本智能体绕开工具发起命令。

## 验证和证据

- `check_symmetry(timeout_s=120)`未找到额外对称操作；没有据此宣称严格证明唯一空间群，也未开展中心对称分支精修。
- `validate_structure(expect_framework=false)`确认单个完整0维分子；其自动“高置信”评分未作为发表质量证据。残峰产生的无序提示未盲目采纳。
- `analyze_packing(blocks=["interactions"], kinds=["hbond"])`提供3条可测量氢键，经EQIV/HTAB进入最终SHELXL CIF。
- `write_outputs`由最终节点产生交付；`run_checkcif`在交付版本3上完成，作业为`.crystalpilot/refine/checkcif/job_20260921_232639_251720_06ee4685/`，原始`model.chk/.ckf/.vrf`及日志保留。
- `analysis/audit_reflections.py`是只读DIALS/ROD诊断脚本，打印目标反射的原始积分、缩放标记和像素统计。工作环境使用给定DIALS Python；输出为UTF-8的`analysis/reflection_audit.txt`。不生成修改后的HKL，不执行索引、积分、缩放或结构精修。

只读核查复现命令，从项目根目录运行：

```powershell
& 'C:\users\<user>\miniforge3\envs\dials\python.exe' -X utf8 'CrystalPilot Results/task_20260921_230205/analysis/audit_reflections.py'
```

引擎自带Python未安装dials，故只读反射诊断改用environment.json给出的DIALS Python。完整dxtbx格式对象初始化曾无诊断退出，重试后仍如此；最终使用同一官方FormatROD的头读取和解压实现，只解码像素并采用已保存探测器面板偏移，核查脚本已成功运行。没有修改工具配置或访问限制。

## 交付边界

本目录保存可直接检查和复算的结构文件、交付来源、验证报告及只读证据。完整大型中间表和计算原始日志仍由本项目`.crystalpilot`审计库持有，未在shell手工拼造另一套结构文件。缺p4p；无溶剂掩膜所以无fab。最终封存不提升为发表级，后续必须处理像素掩膜、重积分、缩放和晶胞误差传播。
