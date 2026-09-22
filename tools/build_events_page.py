"""Build SUBMISSION/运行事件记录.html: an offline page with the inputs, the six conditions, the automatic evaluation
table, the human review verdicts, one screenshot of a CrystalPilot run, the geometry checks of the delivered CIFs
and, for each of the twelve runs, the event stream (tool calls and agent messages), the agent's key statements,
the delivered structure projected in 2D and the verifier's result. Images are embedded as base64; no network."""
import base64
import html
import io
import json
import math
import os
import re
from pathlib import Path

from build_comparison import COND_ZH, RUNS, STATUS_ZH, fmt, pct

ROOT = Path(r"H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921")
CTRL = ROOT / "control"
SUB = CTRL / "SUBMISSION"
OUT = SUB / "运行事件记录.html"
SCREENSHOT = SUB / "images" / "NU1000-CrystalPilot运行中截图.png"
GATEWAY_RE = re.compile(r"https?://[^\s\"']+")
CASE_ZH = {"nu1000": "NU-1000", "alanine": "丙氨酸"}

COND_DESC = {
    "C0N": ("codex exec 0.155.0，无 MCP", "只给一个含 numpy、scipy、matplotlib、pandas、pillow 的 Python 3.12 与 Windows 自带命令；消息前置禁令：不使用任何晶体学软件或程序库，不下载安装，一切算法自行编写。"),
    "P1N": ("claude 2.1.261，无 MCP", "同上。"),
    "C0": ("codex exec 0.155.0，无 MCP", "environment.json 列出 CrystalPilot 自带的 SHELXT、SHELXL、PLATON、Olex2、DIALS 3.30 的路径与用法提示；不给工具、状态、知识、界面。"),
    "P1": ("claude 2.1.261，无 MCP", "同上。"),
    "H0": ("CrystalPilot 工作台，knowledge_mode = tools_only", "73 个专业工具、节点库、工具内校验、交付流程；没有 skill 工具和知识卡片。"),
    "H1": ("CrystalPilot 工作台，knowledge_mode = full", "比纯工具多四个 skill 工具和去标识后的 26 张知识卡片。"),
}

EXPERT = [
    ("纯 Codex、纯 Claude Code", "交付的结构不可用。四个 CIF 都存在结构上的错误，这些错误在 R 值和 checkCIF 警报上没有表现（丙氨酸两轮 R1 为 0.033，键长在常见范围内），只有打开结构逐项检视才能发现。数据处理是临时编写的近似算法，没有任何经过验证的晶体学程序参与；交付的 CIF 缺少存档所需的实验与精修信息。"),
    ("Codex + CrystalPilot 环境、Claude Code + CrystalPilot 环境", "结构上没有大错。但 Agent 为了得到更好的指标，对原始数据做了不利于科学严谨的处理：丙氨酸一轮（r08）在发现低角强反射落在探测器阴影中后，自行定义探测器掩膜并重新积分，把阴影区从数据记录中抹去，没有把这一仪器问题报告出来留待处理；NU-1000 一轮（r04）的 CIF 缺少掩膜信息，另一轮（r07）没有完成 checkCIF。"),
    ("CrystalPilot 纯工具、CrystalPilot 全功能", "结构解析充分，数据处理没有为了指标而调整。丙氨酸两轮都查到了同一处探测器阴影，都保留原始数据与像素证据、拒绝删点或掩盖，并把结果按诊断交付封存；NU-1000 两轮在掩膜收敛性、电子数是否可信、几何是否需要约束上的判断与人类专家一致。数据解析程度与真实性达到与人类专家一致的水平，在数据处理的诚实上更优。全过程在界面上逐节点可见、可引用。"),
]

KEY = {
    "r10_C0N_nu1000": [("最后说明", "已恢复完整三维框架并保存可复现结果，但尚未达到临近发表级精修。全数据精修 R1（I > 2σ）= 0.1965；低角残差、Zr 散射校准和节点质子归属仍未解决，未建立孔道客体模型。"),
                       ("README_zh.md", "从实测强度计算 Patterson 图；检验由重原子矢量提出的 Zr 特殊位置假设。1.0 至 3.0 Å 高角对照（4193 个反射）R1 0.123。")],
    "r09_P1N_nu1000": [("消息 9/23", "训练反射变好，但未参与拟合的反射明显变差，因此这组校正不会进入最终结构。"),
                       ("最后说明", "晶胞内 3 个 Zr₆ 节点、6 个 C₄₄O₈ 配体，形成三维网络，所有 Zr 八配位。"),
                       ("README_ZH.md", "所有晶体学读取、对称运算、合并、Patterson/Fourier 求解、精修及验证算法均在本轮编写。")],
    "r04_C0_nu1000": [("最后说明", "孔道客体未作分子指认，仅以独立散射掩膜处理。"), ("checkCIF", "PLAT609：CIF 缺少 SQUEEZE/MASK 说明；去掩膜 R1 0.1817。")],
    "r07_P1_nu1000": [("最后说明", "完整 checkCIF 验证尚未完成。"), ("复算", "主模型分辩率 1.1 Å（其他五轮 1.0 Å），wR2 0.316。")],
    "r05_H0_nu1000": [("seq 1384，22:36", "掩膜后的初步精修将 R1 降至约 0.0845，但更新掩膜时电子数明显变化且未收敛，因此这还不能作为定稿结果。孔道电子数目前不能用于推算溶剂含量。"),
                      ("seq 1479", "无约束模型中的苯环键长出现约 1.36 至 1.50 Å 的不均衡。我会在独立分支加入芳香环等距、平面性和配体刚性键位移约束，并与无约束结果比较；不改变原子数或占有率。")],
    "r01_H1_nu1000": [("seq 1717，19:43", "删除回峰检验表明，孔道中的不少散射是真实的，但现有“原子”模型出现严重短接触、漂移和极大的 ADP，不能作为化学结构保留。我将保留这套诊断记录，另建纯框架分支。"),
                      ("seq 2040", "末端氧的两个分裂试验均被精修否决，已回到单个位点模型。孔道散射掩膜收敛后，同条件 R1 从 0.1796 降到 0.0845；这支持保留散射校正，但不支持给孔道内容指定分子式。")],
    "r11_C0N_alanine": [("最后说明", "已用指定 Python 自编流程完成解析与精修，并通过从 616 张原始帧开始的完整复现。"), ("METHODS.md", "solve.py 使用固定种子、对称约束的电荷翻转；差值 Fourier 只用实测反射。")],
    "r12_P1N_alanine": [("REPORT_zh.md", "早期电荷翻转未成功；首次完整重跑还暴露了 C/N 指认混淆，已通过至少两个搜索起点和显式指认比较解决。"),
                        ("METHODS_zh.md", "有效求解使用 search_molecule.py：通用丙氨酸连接关系、理想键几何、3 平移 + 3 旋转 + 1 羧基扭转角的差分进化。")],
    "r03_C0_alanine": [("消息 21/28", "PLATON 的当前 Windows 版本停留在图形对话框，尚未产生验证报告。"), ("消息 23/28", "PLATON 因命令行解析错误未完成验证，该限制会保留在报告中。")],
    "r08_P1_alanine": [("消息 17/29", "若干理论上很强的反射被积分为几乎零强度，且集中在某些低角度方向。我会先回查这些反射在原始帧上的位置和像素。"),
                       ("消息 18/29", "原始像素检查已证实存在遮挡。我将根据探测器坐标加入保守的遮挡掩膜后重新积分。"),
                       ("人工审核", "阴影区是采集环节的仪器问题；由 Agent 自行掩膜后重积分让指标变好，却把原始问题从数据记录中移除。")],
    "r06_H0_alanine": [("seq 2176，23:18", "原始帧核查已确认：一条异常 022 观测的预测位置附近，11×11 像素区域全部为零，支持探测器遮挡或屏蔽区污染的判断。当前工具未提供像素掩膜重积分入口，我会将这一限制与误剔除的强观测一起列入交付问题，不用删点掩盖它。")],
    "r02_H1_alanine": [("seq 1960，20:21", "同一等价反射有约 9500 计数的清晰斑点，却被缩放判为离群；接近零、误差很小的测量反而保留下来。这可能与遮挡区域未屏蔽有关。"),
                       ("seq 2042", "第二条还原路线在打开实验时超时，未生成可用结果；现有帧工具也没有暴露遮挡掩膜控制。我会保留这项限制，但不会把低 R 值当作已解决数据异常的证据。"),
                       ("seq 2572，20:29", "原始采集设置明确启用了喷嘴和晶体支架阴影处理，而本次 DIALS 记录只有尺度及衰减校正。我会保存原始帧的对应像素证据，并在最终文件中明确标注：候选结构已收敛，但尚不具备发表级数据验证。")],
}

COLORS = {"C": "#404040", "N": "#3050f8", "O": "#ff0d0d", "H": "#c8c8c8", "Zr": "#3fbfbf", "Zn": "#7d80b0", "Cl": "#1ff01f", "S": "#e6d200", "F": "#90e050", "Br": "#a62929", "Cu": "#c88033"}
SIZES = {"H": 6, "C": 22, "N": 24, "O": 24, "Zr": 60}


def b64(path: Path, mime: str) -> str:
    return f"data:{mime};base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def esc(s) -> str:
    return html.escape(str(s if s is not None else ""))


def cut(s, n: int) -> str:
    """Redact first, then truncate, so a cut never splits a token the redaction must match."""
    return redact(str(s or ""))[:n]


def redact(s: str) -> str:
    s = s.replace(ROOT.name, "scxrd-agent-eval-20260921")
    s = GATEWAY_RE.sub("<url>", s)
    s = re.sub(r"cpl-[0-9a-f]{8,}", "***", s)
    user = os.environ.get("USERNAME", "")
    return re.sub(r"(?i)users([\\/]+)" + re.escape(user), r"users\1<user>", s) if user else s


def structure_png(cif: Path) -> str | None:
    """Unit-cell content of the delivered CIF (all symmetry images of the independent atoms), projected along c
    and along a, with the cell outline."""
    import gemmi
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    st = gemmi.read_small_structure(str(cif))
    if not st.sites:
        return None
    cell = st.cell
    sg = gemmi.find_spacegroup_by_name(st.spacegroup_hm) if st.spacegroup_hm else None
    ops = sg.operations() if sg else gemmi.GroupOps([gemmi.Op("x,y,z")])
    seen = set()
    atoms = []
    for s in st.sites:
        for op in ops:
            fx, fy, fz = op.apply_to_xyz([s.fract.x, s.fract.y, s.fract.z])
            key = (s.element.name, round(fx % 1.0, 3), round(fy % 1.0, 3), round(fz % 1.0, 3))
            if key in seen:
                continue
            seen.add(key)
            p = cell.orthogonalize(gemmi.Fractional(key[1], key[2], key[3]))
            atoms.append((s.element.name, p.x, p.y, p.z))
    order = {"H": 0, "C": 1, "N": 2, "O": 3}
    atoms.sort(key=lambda t: order.get(t[0], 9))
    corners = [cell.orthogonalize(gemmi.Fractional(x, y, z)) for x, y, z in ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.4))
    for ax, (i, j), title, cell_idx in ((axes[0], (0, 1), "沿 c 轴投影（一个晶胞）", (0, 1, 2, 3)), (axes[1], (1, 2), "沿 a 轴投影（一个晶胞）", (0, 3, 7, 4))):
        xs = [a[1 + i] for a in atoms]
        ys = [a[1 + j] for a in atoms]
        cs = [COLORS.get(a[0], "#ff69b4") for a in atoms]
        ss = [SIZES.get(a[0], 26) for a in atoms]
        poly = [(getattr(corners[k], "xyz"[i]), getattr(corners[k], "xyz"[j])) for k in cell_idx] + [(getattr(corners[cell_idx[0]], "xyz"[i]), getattr(corners[cell_idx[0]], "xyz"[j]))]
        ax.plot([p[0] for p in poly], [p[1] for p in poly], color="#888", lw=0.8)
        ax.scatter(xs, ys, c=cs, s=ss, edgecolors="k", linewidths=0.25)
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Å")
        ax.set_ylabel("Å")
    els = sorted({a[0] for a in atoms})
    fig.suptitle(f"{cif.name}：{st.spacegroup_hm}，晶胞内 {len(atoms)} 个原子，元素 {', '.join(els)}", fontsize=10)
    buf = io.BytesIO()
    fig.tight_layout()
    fig.savefig(buf, format="png", dpi=110)
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")


def timeline(cond: str, evid: Path) -> list:
    """(elapsed_min or None, label, detail) rows of the science phase from the native event stream."""
    rows = []
    if cond in ("H0", "H1"):
        p = evid / "raw_events" / "workbench_events.jsonl"
        t_start = None
        for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
            try:
                rec = json.loads(line)
                e = rec["event"]
            except Exception:
                continue
            k = e.get("kind")
            ts = rec.get("ts")
            if k == "user_message":
                t_start = ts if t_start is None else t_start
                if "inputs/" in (e.get("text") or ""):
                    t_start = ts
            if t_start is None:
                continue
            if k in ("tool_started", "command_started"):
                lab = e.get("tool") or e.get("name") or ("shell" if k == "command_started" else "tool")
                det = cut(json.dumps(e.get("args") or e.get("params") or e.get("detail") or {}, ensure_ascii=False), 160)
                rows.append(((ts - t_start) / 60, lab, det))
            elif k == "agent_message":
                rows.append(((ts - t_start) / 60, "message", cut((e.get("text") or ""), 300)))
            elif k in ("approval_request", "background_turn_completed", "turn_completed", "turn_failed", "compaction_started"):
                rows.append(((ts - t_start) / 60, k, ""))
    elif cond in ("C0", "C0N"):
        p = evid / "raw_events" / "codex_events.jsonl"
        for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
            try:
                e = json.loads(line)
            except Exception:
                continue
            it = e.get("item") or {}
            if e.get("type") == "item.completed":
                if it.get("type") == "command_execution":
                    rows.append((None, "shell", cut(it.get("command"), 200) + f" -> exit {it.get('exit_code')}"))
                elif it.get("type") == "agent_message":
                    rows.append((None, "message", cut((it.get("text") or ""), 300)))
                elif it.get("type") == "reasoning":
                    rows.append((None, "reasoning", cut((it.get("text") or ""), 160)))
                elif it.get("type") == "file_change":
                    rows.append((None, "file_change", cut(json.dumps(it.get("changes") or it, ensure_ascii=False), 160)))
    else:
        p = evid / "raw_events" / "claude_stream.jsonl"
        for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
            try:
                e = json.loads(line)
            except Exception:
                continue
            if e.get("type") == "assistant":
                for blk in (e.get("message") or {}).get("content") or []:
                    if not isinstance(blk, dict):
                        continue
                    if blk.get("type") == "tool_use":
                        rows.append((None, blk.get("name"), cut(json.dumps(blk.get("input") or {}, ensure_ascii=False), 200)))
                    elif blk.get("type") == "text" and blk.get("text"):
                        rows.append((None, "message", cut(blk["text"], 300)))
    return rows


def main() -> int:
    summary = json.loads((CTRL / "analysis" / "summary.json").read_text(encoding="utf-8"))
    rows = {r["rid"]: r for r in json.loads((CTRL / "analysis" / "comparison_rows.json").read_text(encoding="utf-8"))}
    checks = json.loads((CTRL / "analysis" / "structure_checks.json").read_text(encoding="utf-8"))
    parts = []
    parts.append(f"""<!doctype html><html lang="zh"><head><meta charset="utf-8"><title>十二轮运行事件记录</title>
<style>body{{font-family:'Microsoft YaHei',Segoe UI,sans-serif;margin:24px auto;padding:0 16px;color:#222;max-width:1240px;line-height:1.55}}h1,h2,h3{{font-weight:600}}
h1{{font-size:24px}}h2{{font-size:19px;margin-top:36px;border-bottom:1px solid #dde3ea;padding-bottom:6px}}h3{{font-size:16px}}
table{{border-collapse:collapse;font-size:13px}}td,th{{border:1px solid #ccc;padding:4px 8px;vertical-align:top}}th{{background:#E7ECF3}}
.run{{border:1px solid #ddd;border-radius:6px;padding:12px 14px;margin:18px 0}}.tl{{font-size:12px;max-height:440px;overflow:auto;border:1px solid #eee;padding:6px}}
.tl div{{padding:2px 0;border-bottom:1px dotted #eee}}.msg{{background:#f7f7ff}}.tag{{display:inline-block;padding:1px 6px;border-radius:3px;background:#eef;margin-right:6px;font-family:Consolas,monospace;font-size:11px}}
img{{max-width:100%}}details summary{{cursor:pointer;color:#246}}.small{{color:#666;font-size:12px}}.quote{{border-left:3px solid #9ab;padding:4px 10px;margin:6px 0;background:#fafcff}}
.shot{{margin:14px 0}}.shot img{{border:1px solid #ddd}}.verdict td:first-child{{white-space:nowrap;font-weight:600}}.c{{text-align:center}}</style></head><body>
<h1>十二轮运行事件记录：通用编程 Agent 与 CrystalPilot 在单晶结构解析上的对照实验</h1>""")
    # 1 inputs
    parts.append("<h2>1　输入</h2><table><tr><th>案例</th><th>输入</th><th>预览</th></tr>")
    ala_jpg = sorted((ROOT / "staging" / "alanine" / "inputs" / "alanine" / "frames" / "jpg").glob("*.jpg"))
    prev = f'<img src="{b64(ala_jpg[0], "image/jpeg")}" style="max-width:260px">' if ala_jpg else ""
    ala_m = json.loads((ROOT / "staging" / "alanine" / "manifest.json").read_text(encoding="utf-8"))
    parts.append(f"<tr><td>丙氨酸</td><td>{len(ala_m)} 个文件（646 帧 .rodhypix，Rigaku HyPix-Arc，Mo Kα），{sum(v['size'] for v in ala_m.values())/1e6:.0f} MB；采集元数据；样品化学式 C3H7NO2。任务提示：用原始衍射数据完成结构解析与精修，得到接近发表质量的结构，保留可复现文件并说明未解决问题。</td><td>{prev}<div class='small'>采集时的晶体照片</div></td></tr>")
    nu_ins = (ROOT / "staging" / "nu1000" / "inputs" / "start.ins").read_text(encoding="utf-8", errors="replace")
    parts.append(f"<tr><td>NU-1000</td><td>start.hkl（337 815 条反射，波长 0.68883 Å）、start.ins（右；只有晶胞、P6/mmm 对称操作、SFAC C H N O 与四个占位原子）、一张合成条件截图（Agent 从中读取金属种类）。范围句：只解析出临近发表级的框架结构，孔道内部的溶剂分子及客体无需处理。</td><td><pre style='font-size:11px;margin:0'>{esc(nu_ins)}</pre></td></tr></table>")
    # 2 conditions
    parts.append("<h2>2　六个条件</h2><p>同一模型 gpt-6-astra（Codex 内核与 CrystalPilot 推理档 xhigh；Claude Code 走网关的 Anthropic 兼容路由，推理档为网关默认值），同一任务提示，同一预算（NU-1000 180 分钟、丙氨酸 120 分钟），无人值守，串行运行。</p><table><tr><th>条件</th><th>内核</th><th>给了什么</th></tr>")
    for cond in ("C0N", "P1N", "C0", "P1", "H0", "H1"):
        k, d = COND_DESC[cond]
        parts.append(f"<tr><td>{esc(COND_ZH[cond])}<br><span class='small'>代码 {cond}</span></td><td>{esc(k)}</td><td>{esc(d)}</td></tr>")
    parts.append("</table>")
    # 3 automatic evaluation
    parts.append("<h2>3　自动评价</h2><p class='small'>R1 独立复算：带 SHELX 文件的轮次为 RES 零周期；纯 Codex/Claude Code 轮次为 CIF 模型转 SHELXL、固定原子只精修标度。checkCIF 为本地 PLATON，A 列括号内为剔除元数据类后的实质警报数。分数区间为预注册评分表的已确认下界至上界，其中两项对纯 Codex/Claude Code 条件结构性不利。</p>")
    parts.append("<table><tr><th>轮</th><th>条件</th><th>案例</th><th>R1 报告</th><th>R1 复算</th><th>差值</th><th>去掩膜 R1</th><th>wR2</th><th>分辩率 Å</th><th>完整度</th><th>精修反射</th><th>checkCIF A（实质）/B/C</th><th>科学状态</th><th>分数区间</th><th>工具调用</th><th>用时 min</th><th>输出 token</th></tr>")
    for rid in RUNS:
        r = rows[rid]
        parts.append(f"<tr><td>{esc(rid)}</td><td>{esc(r['cond_zh'])}</td><td>{CASE_ZH[r['case']]}</td><td class='c'>{fmt(r['r1_rep'])}</td><td class='c'>{fmt(r['r1_rec'])}</td><td class='c'>{fmt(r['r1_diff'])}</td><td class='c'>{fmt(r['nomask'])}</td><td class='c'>{fmt(r['wr2'], 3)}</td><td class='c'>{fmt(r['dmin'], 2)}</td><td class='c'>{pct(r['compl'])}</td><td class='c'>{fmt(r['nrefl'], 0)}</td><td class='c'>{r['A']}（{r['A_subst']}）/{r['B']}/{r['C']}</td><td>{esc(STATUS_ZH.get(r['status'], r['status']))}</td><td class='c'>{fmt(r['lo'], 1)} 至 {fmt(r['hi'], 1)}</td><td class='c'>{r['tools']}</td><td class='c'>{r['wall_min']}</td><td class='c'>{r['out_tok']:,}</td></tr>")
    parts.append("</table>")
    # 4 expert review
    parts.append("<h2>4　人工结构审核</h2><p>审核对象是十二轮交付的主 CIF 及其反射文件、说明文件和过程记录；审核人为晶体学研究者；审核在自动评价完成后进行，逐个打开结构检查。有一类结构错误在 R 值和 checkCIF 警报上不表现，只能由人看结构发现，所以这一层与第 3 节的自动指标并列。</p><table class='verdict'><tr><th>条件层</th><th>结论</th></tr>")
    for a, b in EXPERT:
        parts.append(f"<tr><td>{esc(a)}</td><td>{esc(b)}</td></tr>")
    parts.append("</table>")
    # 5 screenshot of a live run
    if SCREENSHOT.exists():
        parts.append("<h2>5　CrystalPilot 运行中</h2><div class='shot'>"
                     f"<img src='{b64(SCREENSHOT, 'image/png')}'>"
                     "<div class='small'>r05（CrystalPilot 纯工具，NU-1000）运行中的工作台：中栏为 Agent 关于掩膜稳定性、约束与键表修正的判断和逐条工具行（溶剂掩膜、SHELXL 精修、空间群核查、反射数据体检、checkCIF），右栏为节点 n0023 的结构卡（R1 0.0853、wR2 0.2563、GooF 0.91、完整度 88.2%）与框架结构视图。</div></div>")
    # 6 geometry checks
    parts.append("<h2>6　交付结构的几何检查</h2><p class='small'>对每个交付的主 CIF 展开对称、按共价半径搜索键，给出各元素对的键长范围、非氢原子 Ueq 范围和是否有非正定位移参数。十二个 CIF 在这些检查上都没有异常；人工审核在纯 Codex/Claude Code 四轮中发现的结构错误不在这些检查的覆盖范围内。</p>")
    parts.append("<table><tr><th>轮</th><th>条件</th><th>空间群</th><th>独立原子（其中 H）</th><th>键长范围 Å（键数）</th><th>非氢 Ueq 范围</th><th>非正定</th></tr>")
    for rid in RUNS:
        c = checks.get(rid) or {}
        bs = c.get("bond_summary") or {}
        bond = "；".join(f"{k} {v['min']:.3f} 至 {v['max']:.3f}（{v['n']}）" for k, v in sorted(bs.items()))
        parts.append(f"<tr><td>{esc(rid)}</td><td>{esc(rows[rid]['cond_zh'])}</td><td>{esc(c.get('space_group'))}</td><td class='c'>{esc(c.get('n_sites'))}（{esc(c.get('n_H'))}）</td><td>{esc(bond)}</td><td class='c'>{esc(c.get('ueq_nonH_min'))} 至 {esc(c.get('ueq_nonH_max'))}</td><td class='c'>{'无' if not c.get('npd') else esc(c.get('npd'))}</td></tr>")
    parts.append("</table>")
    # 7 per run
    parts.append("<h2>7　逐轮事件记录</h2><p class='small'>每轮给出 Agent 的关键原话（工作台事件流的 seq 或原生事件流的消息序号）、交付结构在一个晶胞内的投影、完整事件流（默认折叠；工作台条件带自任务开始的分钟数）、最后说明与独立复算结果。</p>")
    for rid in RUNS:
        r = rows[rid]
        s = summary["runs"][rid]
        row = s["row"]
        evid = CTRL / "runs" / rid
        parts.append(f"<div class='run'><h3>{esc(rid)}　{esc(r['cond_zh'])} / {CASE_ZH[r['case']]}</h3>")
        parts.append(f"<p><span class='tag'>科学状态 {esc(STATUS_ZH.get(r['status'], r['status']))}</span><span class='tag'>R1 报告 {fmt(r['r1_rep'])} / 复算 {fmt(r['r1_rec'])}</span><span class='tag'>用时 {r['wall_min']} min</span><span class='tag'>工具调用 {r['tools']}</span><span class='tag'>退出 {esc(row.get('exit_reason'))}</span><span class='tag'>自动干预 {esc(row.get('automatic_interventions'))}</span><span class='tag'>人工干预 0</span></p>")
        for pos, q in KEY.get(rid, []):
            parts.append(f"<div class='quote'><span class='small'>{esc(pos)}</span><br>{esc(q)}</div>")
        if s.get("fatal"):
            parts.append("<p><b>复算/校验发现：</b>" + "; ".join(esc(x) for x in s["fatal"]) + "</p>")
        cif_rel = ((s.get("delivery") or {}).get("picked") or {}).get(".cif")
        if cif_rel:
            arena = Path(json.loads((evid / "manifest.json").read_text(encoding="utf-8"))["arena"])
            try:
                img = structure_png(arena / cif_rel)
                if img:
                    parts.append(f"<img src='{img}'>")
            except Exception as e:  # noqa: BLE001
                parts.append(f"<p class='small'>结构投影不可用：{esc(type(e).__name__)}: {esc(e)}</p>")
        try:
            tl = timeline(r["cond"], evid)
        except Exception as e:  # noqa: BLE001
            tl = [(None, "error", f"timeline unavailable: {type(e).__name__}")]
        parts.append(f"<details><summary>事件流（{len(tl)} 条）</summary><div class='tl'>")
        for t, lab, det in tl[:600]:
            tt = f"{t:6.1f} min" if isinstance(t, (int, float)) else ""
            cls = " class='msg'" if lab == "message" else ""
            parts.append(f"<div{cls}><span class='tag'>{esc(tt)} {esc(lab)}</span>{esc(redact(det))}</div>")
        parts.append("</div></details>")
        if s.get("final_text_tail"):
            parts.append(f"<details><summary>Agent 的最后说明</summary><pre style='white-space:pre-wrap;font-size:12px'>{esc(redact(s['final_text_tail']))}</pre></details>")
        v = s.get("verdict") or {}
        parts.append(f"<details><summary>独立复算与评分</summary><pre style='font-size:12px'>{esc(json.dumps({k: v.get(k) for k in ('r1_reported', 'r1_recomputed', 'r1_abs_difference', 'r1_definition', 'fatal_findings', 'score_lower_bound', 'score_upper_bound', 'unknown_items')}, ensure_ascii=False, indent=1))}</pre></details>")
        parts.append("</div>")
    parts.append("</body></html>")
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print("written", OUT, OUT.stat().st_size, "bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
