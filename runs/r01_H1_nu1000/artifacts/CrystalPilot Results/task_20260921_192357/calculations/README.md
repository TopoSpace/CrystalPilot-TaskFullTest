# 计算记录与复算说明

本目录只包含本轮当前项目的记录；原始inputs未修改。主交付在上一级。所有求解、改模、精修及结构导出通过CrystalPilot工具完成；shell仅用于读取、分析、复制和报告整理。

## 目录

- `shelxt/`：作业job_20260921_192553的完整输入/输出，含job.ins/hkl、job.lxt、stdout、各候选res/hkl。job_a为采纳的P6/mmm解；组成Zr6 C88 O32、z=6只是相位辅助猜测，不是最终化学式/Z。
- `shelxl/<job-id>/`：本轮全部30个SHELXL作业的.ins/.res/.lst/.fcf/.fab（实际存在者），包括失败/拒绝作业。未重复存每份相同的job.hkl，使用上一级`final.hkl`的原始观测。零字节res代表失败结果，不能当结构读取。
- `nodes/`：n0000至n0059的工具节点快照，含model.res/cif、node.json、峰表及存在时的掩膜快照。节点CIF是工作格式，不替代主交付带esd/实验字段的final.cif。
- `checkcif/`：最终交付版本3的完整PLATON作业及报告。
- `source_inputs/`：原始start.ins和合成图的只读来源副本；原始HKL由主交付final.hkl保留。
- `investigation.json`、`ghost_ledger.json`、`trial_ledger.json`：本轮目标、已否决配置和假设检验账本。

## 关键作业

|作业后缀（日期均20260921）|用途/节点|
|---|---|
|194003、194044|初始各向同性/各向异性框架及未解释孔道位点|
|194058至194224|ghost_test参考与12组删除/回峰计算|
|194457|无掩膜纯框架基线n0037|
|194608|无掩膜末端O分裂n0041，不采纳|
|194717|首次收敛掩膜后n0043，同权重对照|
|195020|刷新但未收敛掩膜n0045，不采纳|
|195137|Zr2分裂不稳定，失败记录|
|195220|局部配体试验过程中输入/氢连接检查记录|
|195250|局部配体分裂n0051，不采纳|
|195334、195334_w1|权重采纳环，最终WGHT0.1444,0|
|195425|Z=3下的收敛模型n0054|
|195520|FREE同名原子指令被拒绝，未改节点|
|195821|修正EQIV/FREE后的最终n0055，主CIF/FCF/INS来源|
|200042|当前掩膜末端O分裂n0059，负占有率，不采纳|
|200134|n0055的SHELXL check，复现最终指标|

若原始`.lst`与文字小数舍入有细微差异，以作业`.lst`及主交付CIF真实字段为准。所有非最终分支都是诊断，不可仅按R较小选成最终结构。

## 人工复算

在工作台内继续时，应从n0055/主交付模型恢复并用CrystalPilot类型化工具精修，保持节点审计。下列命令只说明离线人工核对方式，本轮没有绕过工具运行结构引擎。

1. 将主交付`final.ins`、`final.hkl`、`final.fab`放在一个工作副本目录，避免覆盖封存文件。
2. 用环境提供的SHELXL-2019/3在该副本目录运行：

```powershell
& 'H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921\toolbox\shelx\shelxl.exe' final
```

预期R1=0.0849、wR2=0.2526、GooF=0.912。`.ins`已含SHEL 999 1.000、正确DISP、WGHT、AFIX、ABIN和EQIV/FREE；缺失final.fab会破坏掩膜复算。最终INS含155参数、0 restraints，不需要手动猜软件默认组成。

复算历史作业时，在相应目录的副本中把主`final.hkl`复制为`job.hkl`，并保留该作业原有的`job.fab`（若INS有ABIN），运行`shelxl job`。无掩膜基线194457不需要FAB。SHELXT候选目录有其自己的输入HKL；不要把较低对称候选输出误当最终结构。

只读分析脚本`../analysis/audit_delivery.py`用gemmi/cctbx解析最终文件，明确HKLF4，不改任何模型，输出缺失反射等价类、负观测数量和原子Ueq表。结果已存为audit_delivery.json。

`final.res`适合继续模型；`final.ins`是生成发表格式CIF/FCF的确切输入。无p4p，因为输入未提供，不能由结构假造。CIF中的未知温度/仪器等也不应在复算时自动当作默认实验事实。
