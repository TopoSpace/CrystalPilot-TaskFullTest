# CrystalPilot-TaskFullTest

通用编程 Agent 与晶体学 Agent 工作台 CrystalPilot 在单晶 X 射线衍射结构解析任务上的对照实验：两组真实数据、六个条件、十二轮无人值守运行的完整证据。

An English summary is at the end of this file.

---

## 1　实验概要

| 项目 | 内容 |
|:---:|---|
| 目的 | 检验通用编程 Agent（Codex、Claude Code）直接执行单晶结构解析的能力，以及 CrystalPilot 提供的工具、状态、纪律与交互带来的变化 |
| 时间 | 2026 年 9 月 21 日 19:24 至 22 日 05:27（北京时间），十二轮串行 |
| 模型 | gpt-6-astra，全部条件同一模型 |
| 案例 | 丙氨酸原始衍射图像（646 帧）；锆基 MOF NU-1000 的反射数据 |
| 条件 | 纯 Codex、纯 Claude Code、Codex + CrystalPilot 环境、Claude Code + CrystalPilot 环境、CrystalPilot 纯工具、CrystalPilot 全功能 |
| 评价 | 独立程序复算 R1、本地 PLATON checkCIF、逐条命令审计、晶体学研究者的人工结构审核 |
| 干预 | 自动干预 0 次，人工科学干预 0 次 |

---

## 2　任务

单晶 X 射线衍射结构解析从一颗晶体的衍射图像出发，得到一个可以发表和存档的化学结构。链条上的每个环节都有参数选择，前一步的错误会在两三步之后以变差的统计量出现。

| 环节 | 做什么 | 输出 |
|:---:|---|---|
| 数据还原 | 找衍射点、指标化、积分、缩放合并 | 反射数据（HKL） |
| 空间群判定 | 由对称性与系统消光确定空间群 | 空间群 |
| 结构求解 | 用直接法、电荷翻转或 Patterson 方法恢复相位 | 初始原子模型 |
| 精修 | 最小二乘调整坐标、位移参数、占有率；处理无序、溶剂、氢原子 | 收敛的模型 |
| 校验 | R1、wR2、GooF、残余电子密度；checkCIF 的 A/B/C 级警报 | 校验报告 |
| 交付 | 模型、反射数据、精修记录写成 CIF | CIF 与配套文件 |

一个结构的质量由 R1（小分子发表水平 0.05 以下，多孔框架 0.10 以下）、wR2、GooF、残峰和 checkCIF 的 A 级警报数衡量；结构是否正确，最终由晶体学研究者打开 CIF 检查决定，因为有一类结构错误不反映在任何指标或警报上。

---

## 3　两个案例

| 案例 | 输入 | 化学与对称性 | 难点 | 预算 |
|:---:|---|---|---|:---:|
| 丙氨酸 | 646 张原始衍射图像（Rigaku HyPix-Arc，Mo Kα）及采集元数据，110 MB；化学式 C₃H₇NO₂ | L-丙氨酸，P2₁2₁2₁，晶胞 5.79 × 5.97 × 12.30 Å，6 个非氢原子 | 从图像开始的完整链条；采集中存在探测器阴影 | 120 分钟 |
| NU-1000 | `start.hkl`（337 815 条反射，波长 0.68883 Å）、`start.ins`（只有晶胞、P6/mmm 对称操作、SFAC C H N O 与四个占位原子）、一张合成条件截图 | Zr₆ 节点 + TBAPy 配体，P6/mmm，a = 39.19 Å，c = 16.61 Å，V = 22 093 Å³ | 大晶胞高对称框架；孔道溶剂散射；弱高角数据；金属种类只在截图中 | 180 分钟 |

NU-1000 的任务范围句："只解析出临近发表级的框架结构，孔道内部的溶剂分子及客体无需处理。"输入文件的公开范围见 `inputs/`：NU-1000 的反射数据与起始文件在仓库内；丙氨酸原始图像作为 Release 附件提供，仓库内保留其哈希清单。

---

## 4　六个条件

| 条件 | 内核 | 给了什么 |
|:---:|---|---|
| 纯 Codex | codex exec 0.155.0，无 MCP | 一个只含 numpy、scipy、matplotlib、pandas、pillow 的 Python 3.12 和 Windows 自带命令；消息前置禁令：不使用任何晶体学软件或程序库，不下载安装，一切算法自行编写 |
| 纯 Claude Code | claude 2.1.261，无 MCP | 同上 |
| Codex + CrystalPilot 环境 | codex exec 0.155.0，无 MCP | environment.json 列出 CrystalPilot 自带的 SHELXT、SHELXL、PLATON、Olex2、DIALS 3.30 的路径与用法提示；不给工具、状态、知识、界面 |
| Claude Code + CrystalPilot 环境 | claude 2.1.261，无 MCP | 同上 |
| CrystalPilot 纯工具 | CrystalPilot 工作台，knowledge_mode = tools_only | 73 个专业工具、节点库、工具内校验、交付流程 |
| CrystalPilot 全功能 | CrystalPilot 工作台，knowledge_mode = full | 比纯工具多四个 skill 工具和去标识后的 26 张知识卡片 |

Codex 内核与 CrystalPilot 的推理档位为 xhigh；Claude Code 通过同一网关的 Anthropic 兼容接口调用同一模型，该接口不接受推理档位参数，推理档位为网关默认值。每个条件在每一例上收到相同的任务材料和一句短任务提示（`freeze/prompts/`），各轮提示的 SHA-256 在 `runs/<run_id>/manifest.json`。

---

## 5　实验流程

| 步骤 | 做法 | 仓库中的记录 |
|:---:|---|---|
| 1 环境冻结 | CrystalPilot 按提交 `327e1ef` 导出为冻结副本并另起服务实例；记录内核、模型、科学软件、机器 | `environment/environment.json` |
| 2 输入封存 | 复制原始数据并记录每个文件的哈希；丙氨酸目录中 206 个携带还原结果或数据库检索结果的文件被排除并逐个记录理由 | `provenance/alanine_staging.json`、`nu1000_staging.json`、`inputs/` |
| 3 知识快照去标识 | 全功能模式的 26 张知识卡片中与 NU-1000 直接相关的三处内容改写，前后哈希留存 | `provenance/knowledge_redactions.json`、`freeze/knowledge_snapshot_manifest.json` |
| 4 计划冻结 | 带软件的八轮在首个结果前冻结，顺序由种子 20260921 生成；纯产品四轮在八轮结果出来后追加并冻结 | `freeze/plan.json`、`freeze/plan_pure.json`、`freeze/protocol_frozen.json` |
| 5 无人值守运行 | 运行器按硬时限启动内核或工作台回合，记录原生事件流、命令、干预与用量；总控只做确定性的环境管理 | `tools/run_trial.py`、`tools/run_all.py`；`runs/<run_id>/` |
| 6 封存 | 每轮结束立即记录参试目录全部文件的哈希与输入是否被改动 | `runs/<run_id>/seal_manifest.json` |
| 7 独立复算 | 同一版本 SHELXL 在交付的反射集合上零周期复算 R1；纯产品四轮把 CIF 模型转成 SHELXL 作业只精修标度 | `tools/verify_run.py`、`tools/recompute_from_cif.py`；`verifier/<run_id>/` |
| 8 checkCIF | 本地 PLATON checkCIF，A 级警报分元数据类与实质类 | `verifier/<run_id>/checkcif/` |
| 9 命令审计 | 纯产品四轮的全部命令与文件访问逐条核对禁令 | `tools/audit_pure.py`；`analysis/pure_audit.json` |
| 10 汇总 | 逐轮结果表、对比表、几何检查、节点计数 | `tools/analyze.py`、`tools/build_comparison.py`、`tools/structure_checks.py`；`analysis/` |
| 11 人工结构审核 | 晶体学研究者逐个打开十二个交付的 CIF 及配套文件检查 | 结论见第 6.2 节 |

运行规则全文见 `PROTOCOL.md`。

---

## 6　结果

### 6.1　自动评价

R1 独立复算：带 SHELX 文件的轮次为 RES 零周期，纯产品轮次为 CIF 模型转 SHELXL 只精修标度。GooF 以接近 1 为优。checkCIF A 列括号内为剔除元数据类后的实质警报数。完整表（含去掩膜 R1、Rint、复算途径、差值、分数区间、token 用量）见 `analysis/comparison_tables.md` 与 `analysis/results.csv`。

| 轮 | 条件 | 案例 | R1 报告 | R1 复算 | wR2 | GooF | 分辩率 Å | 完整度 | 精修反射 | checkCIF A（实质）/B/C | 工具调用 | 用时 min |
|---|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| r10 | 纯 Codex | NU-1000 | 0.1965 | 0.1962 | 0.504 | 1.742 | 全数据 | 未记录 | 11289 | 22（16）/7/12 | 63 | 39.7 |
| r09 | 纯 Claude Code | NU-1000 | 0.1750 | 0.1785 | 0.463 | 2.513 | 1.00 | 未记录 | 4395 | 21（15）/5/12 | 159 | 63.0 |
| r04 | Codex + CrystalPilot 环境 | NU-1000 | 0.0827 | 0.0827 | 0.272 | 0.987 | 1.00 | 99.9% | 4395 | 10（3）/2/23 | 174 | 47.5 |
| r07 | Claude Code + CrystalPilot 环境 | NU-1000 | 0.0876 | 0.0876 | 0.316 | 1.021 | 1.10 | 99.9% | 3345 | 9（2）/2/21 | 178 | 45.3 |
| r05 | CrystalPilot 纯工具 | NU-1000 | 0.0853 | 0.0853 | 0.256 | 0.913 | 1.00 | 99.9% | 4395 | 9（2）/1/23 | 93 | 46.5 |
| r01 | CrystalPilot 全功能 | NU-1000 | 0.0849 | 0.0849 | 0.253 | 0.912 | 1.00 | 99.9% | 4395 | 9（2）/1/20 | 140 | 47.4 |
| r11 | 纯 Codex | 丙氨酸 | 0.0326 | 0.0326 | 0.075 | 1.311 | 0.75 | 97.0% | 623 | 3（2）/3/9 | 96 | 54.4 |
| r12 | 纯 Claude Code | 丙氨酸 | 0.0334 | 0.0335 | 0.076 | 1.329 | 0.71 | 86.3% | 657 | 7（3）/3/3 | 199 | 48.4 |
| r03 | Codex + CrystalPilot 环境 | 丙氨酸 | 0.0304 | 0.0304 | 0.073 | 1.075 | 0.76 | 99.2% | 984 | 1（0）/0/4 | 187 | 42.3 |
| r08 | Claude Code + CrystalPilot 环境 | 丙氨酸 | 0.0302 | 0.0302 | 0.073 | 1.095 | 0.77 | 99.8% | 964 | 3（0）/0/3 | 200 | 43.8 |
| r06 | CrystalPilot 纯工具 | 丙氨酸 | 0.0457 | 0.0452 | 0.205 | 1.271 | 0.71 | 87.9% | 1055 | 4（1）/3/20 | 74 | 25.7 |
| r02 | CrystalPilot 全功能 | 丙氨酸 | 0.0442 | 0.0437 | 0.200 | 1.222 | 0.71 | 87.9% | 1055 | 4（1）/3/19 | 111 | 28.8 |

十二轮的 R1 都由独立程序复算：带 SHELX 文件的八轮差值不超过 0.0005，纯产品四轮差值 0.0000 至 0.0035。纯产品四轮的命令审计无违规。十二轮的输入文件都没有被改动。

### 6.2　人工结构审核

审核对象是十二轮交付的主 CIF 及其反射文件、说明文件和过程记录，审核人为本实验的晶体学研究者，审核在自动评价完成后进行。

| 条件层 | 审核结论 |
|:---:|---|
| 纯 Codex、纯 Claude Code | 交付的结构不可用。四个 CIF 都存在结构上的错误，这些错误在 R 值和 checkCIF 警报上没有表现（丙氨酸两轮 R1 为 0.033，键长在常见范围内），只有打开结构逐项检视才能发现。数据处理是临时编写的近似算法，没有任何经过验证的晶体学程序参与；CIF 缺少存档所需的实验与精修信息。 |
| Codex + CrystalPilot 环境、Claude Code + CrystalPilot 环境 | 结构上没有大错。但 Agent 为了得到更好的指标，对原始数据做了不利于科学严谨的处理：丙氨酸一轮（r08）在发现低角强反射落在探测器阴影中后，自行定义探测器掩膜并重新积分，把阴影区从数据记录中抹去，而不是把这一仪器问题报告出来留待处理；NU-1000 一轮（r04）的 CIF 缺少掩膜信息，另一轮（r07）没有完成 checkCIF。 |
| CrystalPilot 纯工具、CrystalPilot 全功能 | 结构解析充分，数据处理没有为了指标而调整。丙氨酸两轮都查到了同一处探测器阴影，都保留原始数据与像素证据、拒绝删点或掩盖，并把结果按诊断交付封存；NU-1000 两轮在掩膜收敛性、电子数是否可信、几何是否需要约束上的判断与人类专家一致。数据解析程度与真实性达到与人类专家一致的水平，在数据处理的诚实上更优；全过程逐节点可见、可引用。 |

### 6.3　三层结果

- **纯产品**：两个原生 Agent 都选择自己编写整套算法（解码探测器格式、指标化、积分、合并、求解、精修、写 CIF，1300 至 2100 行代码），得到的 R1 被独立程序证实，但结构经人工审核不可用。这说明这项任务需要晶体学工作流与工作环境。
- **环境组**：有了 SHELX、PLATON、DIALS，结构不再有大错，指标上是十二轮里最好的两轮（丙氨酸 R1 0.030，完整度 99% 以上）；这些数字来自对原始数据记录的处理。
- **CrystalPilot**：在同一处探测器阴影前记录证据、拒绝改动数据、按诊断交付（R1 0.044 至 0.045，完整度 88%）；NU-1000 上的掩膜与约束判断与专家一致；平均用时 36 至 38 分钟、84 至 126 次工具调用，快于其余条件。

---

## 7　目录结构

```
README.md                 本文
PROTOCOL.md               运行规程全文
LICENSE-CODE              tools/ 的许可证（MIT）
LICENSE-DATA.md           文档、证据与数据的许可证（CC BY 4.0）
checksums.sha256          仓库全部文件的 SHA-256
inputs/
  nu1000/                 start.hkl、start.ins、三个输入文件的哈希
  alanine/                原始图像目录的哈希清单；图像本体见 Release 附件
freeze/
  plan.json               带软件八轮的冻结计划
  plan_pure.json          纯产品四轮的冻结计划（含禁令原文与提前停止规则）
  protocol_frozen.json    冻结时刻各文件的哈希
  prompts/                任务提示与共同运行约定
  tool_names_*.json       两种知识模式的工具名清单
  knowledge_snapshot_manifest.json
environment/              版本、模型、机器、科学软件；纯产品条件的 Python 环境
provenance/               输入封存与排除清单；知识快照去标识记录
runs/<run_id>/
  prompt.txt              该轮收到的完整消息
  final_message.txt       Agent 的最后说明
  interventions.jsonl     自动干预记录
  usage.json              工具调用与 token 用量
  manifest.json           条件、模型、用时、退出原因、提示哈希
  seal_manifest.json      封存时参试目录全部文件的哈希；inputs_tampered
  raw_events/             原生事件流（工作台事件、codex 事件或 claude 流）
  commands/               命令记录
  effective_instructions/ 实际生效的指令文件（工作台条件）
  native_state/           节点索引、知识写入、MCP 日志（工作台条件）
  artifacts/              交付的结构文件、说明、节点库与自编代码（文本类型）
verifier/<run_id>/
  report.json             复算与评分报告
  recompute/              SHELXL 零周期复算作业与输出
  recompute_from_cif*/    CIF 转 SHELXL 复算（纯产品四轮）
  checkcif/               PLATON checkCIF 输出
analysis/
  results.csv             逐轮结果表
  run_index.csv           每轮的证据路径与会话 ID
  comparison_tables.md    自动评价全表
  comparison_rows.json    表格数据
  structure_checks.json   交付 CIF 的键长与位移参数检查
  pure_audit.json         纯产品四轮的命令审计
  node_counts.json        工作台条件的节点数
  summary.json            汇总
tools/                    运行器、总控、复算、审计、汇总脚本
```

运行编号中的条件代码：C0N 纯 Codex、P1N 纯 Claude Code、C0 Codex + 环境、P1 Claude Code + 环境、H0 CrystalPilot 纯工具、H1 CrystalPilot 全功能。

---

## 8　核验与复跑

**核验文件完整性**

```bash
sha256sum -c checksums.sha256
```

**核验一轮的封存与复算**：`runs/<run_id>/seal_manifest.json` 记录了参试目录在封存时刻的全部文件哈希与 `inputs_tampered`（十二轮均为空）；`verifier/<run_id>/report.json` 记录了复算所用文件的哈希、SHELXL 与 PLATON 输出的位置。`verifier/<run_id>/recompute/` 中的 `.ins` 与交付的 `.hkl` 可用同版本 SHELXL 直接重跑。

**复跑一轮**需要：Windows；CrystalPilot 提交 `327e1ef` 的检出（含 `.venv` 与 `vendor`）；DIALS 3.30 conda 环境；SHELXL 2019/3、SHELXT 2018/2、PLATON；Codex 内核 0.155.0（随 CrystalPilot 提供）；Claude Code 2.1.261；一个 OpenAI Responses 兼容的模型网关及其凭据；纯产品条件另需按 `environment/pure_baseline_environment.json` 建立的 Python 3.12 虚拟环境。第三方程序不随本仓库分发。

1. 按 `tools/run_trial.py` 顶部的路径说明建立实验根目录，配置引擎副本、三套隔离的内核目录与凭据钩子。
2. `python tools/stage_inputs.py` 封存输入；`python tools/freeze_plan.py --include-p1 yes` 冻结计划。
3. `python tools/run_trial.py --run-id r01_H1_nu1000` 运行一轮；纯产品轮次加 `--plan freeze/plan_pure.json`。
4. `python tools/verify_run.py --run-id <run_id>` 独立复算；纯产品轮次再运行 `python tools/recompute_from_cif.py --run-id <run_id>`。
5. `python tools/analyze.py`、`python tools/build_comparison.py`、`python tools/audit_pure.py`、`python tools/structure_checks.py` 汇总。

---

## 9　与 CrystalPilot 的关系

CrystalPilot 是 TopoSpace 开发的单晶结构解析 Agent 工作台：以 Codex app-server 为内核，通过 MCP 提供类型化的晶体学工具，把每步计算写入可分支回溯的节点库，交付时附带校验文件、逐条警报解释和说明。本实验使用其私有仓库提交 `327e1ef` 的冻结副本，产品本体不在本仓库中；公开预览见 https://github.com/TopoSpace/CrystalPilot-Preview 。本仓库中的文本已脱敏：模型网关地址、凭据与本机用户名已替换为占位符。

---

## 10　许可证

`tools/` 下的脚本采用 MIT 许可证（`LICENSE-CODE`）；文档、证据记录与数据采用 CC BY 4.0（`LICENSE-DATA.md`）。

---

## English summary

This repository holds the complete evidence of a controlled experiment on single-crystal X-ray structure determination, run on 21 and 22 September 2026. One model (gpt-6-astra) processed two real data sets, a raw-frame L-alanine collection (646 Rigaku HyPix frames) and the reflection data of the Zr-MOF NU-1000, under six conditions: bare Codex and bare Claude Code with no crystallographic software at all; Codex and Claude Code given the software environment assembled by CrystalPilot (SHELXT, SHELXL, PLATON, Olex2, DIALS); and the CrystalPilot workbench in its tools-only and full modes. Twelve unattended runs were sealed, re-computed with SHELXL by an independent script, checked with PLATON checkCIF, audited command by command, and finally inspected structure by structure by a crystallographer. The bare agents wrote their own decoders, integration and refinement code and reached plausible R factors, but the delivered structures contain errors that neither R1 nor checkCIF reveals and are unusable. Given the environment, the structures no longer had major errors, but the agents manipulated the data record for better statistics, masking detector-shadowed regions and re-integrating. CrystalPilot diagnosed the same shadow, refused to alter the data, delivered the result as diagnostic, and exposed every node, metric trajectory and alert for the human to cite and steer. Directory layout, verification and re-run instructions are given above; the protocol is in PROTOCOL.md.
