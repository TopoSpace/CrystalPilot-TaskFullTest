# 计算记录与重现

## 范围

输入只来自当前项目 `inputs/alanine/`。未读取其他运行、历史结构、参考答案或评分目录，未联网检索该样品。全部模型创建/改动、精修与 CIF/HKL/FCF 交付均经 `mcp__crystalpilot` 工具；两份 Python 脚本仅作只读分析。原始 inputs 未修改。

## 工具调用顺序

```text
get_project_brief()
set_experiment()
import_frames(frames_dir="inputs/alanine/frames",
              dials_env="C:/users/<user>/miniforge3/envs/dials",
              min_run_frames=12)
find_spots(nproc=1)
index_frames()
integrate_frames(nproc=1)
scale_and_export(composition="C3 H7 N O2")
create_start_model(composition="C3 H7 N O2", symmetry="auto")
screen_space_groups(laue_group="all", merge_stats=true)
check_symmetry()
run_shelxt(composition="C3 H7 N O2", z=4, adopt=false,
           all_space_groups=true, timeout_s=120,
           search_grace_s=120, phasing_grace_s=120)
assemble_asu()
edit_atoms(reassign original C1 to N)
rename_atoms(map={C1:N1, C3:C1, C4:C3})
set_adp(all six non-H atoms, anisotropic)
run_shelxl(mode="adopt", l_s=8)
add_hydrogens(elements=["C"],
              force_kind={N1:CH3_rotating, C3:CH3},
              bond_lengths={CH3_rotating:0.91, CH3:0.98,
                            tertiary_CH:1.00})
run_shelxl(mode="adopt", l_s=10)
audit_reflection_data()
estimate_resolution()
```

定群工具在会话尚未建立时的首次调用被正常拒绝，自动起模后补做完整筛查；粗解空间群仅作为候选，之后由消光、Laue 合并、SHELXT及结构对称检查交叉检验。没有静默换群。

O1 分支从 n0007 开始，`model_disorder` 将次位置置于残差峰，起始主占有率0.95；SADI/SIMU后10循环得到 n0010，再撤销为 n0011。O2 另从 n0007开始，同样处理，15+30循环得到 n0015，再撤销为 n0016。实际 restraints 与判据见 SUMMARY。

交付分支重新从 n0007 开始。`run_shelxl(mode="adopt_wght", l_s=10, wght_rounds=4)` 得 n0017，实际采用 `WGHT 0.1353 0.0753`。`assemble_asu`确认无需再移动，`check_symmetry`未发现漏对称。`run_shelxl(mode="check", l_s=8)` 复现 R1=0.0442。加入三条 EQIV/HTAB 氢键输出卡后 `run_shelxl(mode="adopt", l_s=8)` 得 n0018，其指标相同。

最终 `write_outputs(status="diagnostic")` 配对真实 SHELXL 作业；在交付版本2的 `final.cif` 上运行本地完整 `run_checkcif(timeout_s=180)`，未外发到 IUCr。所有工具预算均在环境总预算内，未修改环境限制。

## 交付目录

- `final.ins` 与 `final.hkl`：与最终 CIF/FCF 配对的可重启 SHELXL 输入。
- `final.res`：节点 n0018 的模型；`final.cif`含精修统计、几何、氢键和嵌入数据；`final.fcf`为实际结构因子。
- `records/frames/`：DIALS各阶段 expt/refl、缩放 MTZ、导出 HKL/INS、参数状态和日志。
- `records/shelxl/`、`records/shelxt/`：所有本轮精修/求解作业输入和输出，不只保留最好一次。
- `records/nodes/`、`records/data/`、`records/runs/`：本轮节点链、数据快照、工具事件。
- `records/checkcif/`：PLATON原始输出；顶层 `checkcif.json`绑定交付版本2。
- `records/views/`：实际查看的结构图片；`analysis/`为原始像素与积分标记的只读诊断。

原始扫描仍在项目 inputs 中，未重复包装616张原图。DIALS expt 文件中的路径指向本项目 `.crystalpilot/frames/frames_normalized/` 的工具生成链接。移动整个项目后需重新导入原始帧；单独移动此结果目录仍可从 final.ins/hkl 复现精修，但不能凭结果目录独立重做原始积分。

## 接手重现

在本项目 CrystalPilot 会话中可 checkout `n0018` 后使用 `run_shelxl(mode="check")`复核，不会移动节点原子。离开工作台时，标准 SHELXL 的等价复算是以同 basename 的 `final.ins`、`final.hkl`运行 `shelxl final`；本轮没有在 shell 中执行此旁路。重现的是当前诊断模型与数据，不是已修复的数据。

原始求解节点登记为 import，导致工具生成的 CIF 中 `_computing_structure_solution`仍为 `?`。已用 set_experiment记录真实求解方法，但写出器按节点谱系仍覆盖此栏；故在此完整披露，不手工拼接 CIF。零晶胞误差也是导出限制，不代表真实零不确定度。
