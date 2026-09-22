"""Write control/analysis/COMPARISON_12_tables.md: the twelve-run comparison tables (r01-r08 plus the pure
baselines r09-r12; the user-requested rerun is excluded), built only from sealed evidence and verifier
reports. The narrative is written separately in COMPARISON_12.md and includes these tables."""
import json
from pathlib import Path

ROOT = Path(r"H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921")
CTRL = ROOT / "control"
RUNS = ["r10_C0N_nu1000", "r09_P1N_nu1000", "r04_C0_nu1000", "r07_P1_nu1000", "r05_H0_nu1000", "r01_H1_nu1000",
        "r11_C0N_alanine", "r12_P1N_alanine", "r03_C0_alanine", "r08_P1_alanine", "r06_H0_alanine", "r02_H1_alanine"]
COND_ZH = {"H1": "CrystalPilot 全功能", "H0": "CrystalPilot 纯工具", "C0": "Codex + CrystalPilot 环境", "P1": "Claude Code + CrystalPilot 环境",
           "C0N": "纯 Codex", "P1N": "纯 Claude Code"}
ORDER = ["C0N", "P1N", "C0", "P1", "H0", "H1"]
STATUS_ZH = {"candidate_checks_passed": "候选检查通过", "diagnostic_delivery": "诊断交付", "partial_result": "部分结果",
             "no_usable_result": "无可用结果", "not_evaluable": "不可评"}


def fmt(v, nd=4):
    if v is None or v == "":
        return "—"
    try:
        f = float(str(v).split("(")[0])
        return f"{f:.{nd}f}" if nd else f"{f:.0f}"
    except ValueError:
        return str(v)


def pct(v):
    if v is None or v == "":
        return "—"
    try:
        return f"{float(str(v).split('(')[0]) * 100:.1f}%"
    except ValueError:
        return str(v)


def main() -> int:
    rows = []
    for rid in RUNS:
        evid = CTRL / "runs" / rid
        man = json.loads((evid / "manifest.json").read_text(encoding="utf-8"))
        us = json.loads((evid / "usage.json").read_text(encoding="utf-8"))
        rep = json.loads((CTRL / "verifier" / rid / "report.json").read_text(encoding="utf-8"))
        v = rep["verdict"]
        cif = (rep["parsed"].get("cif") or {}).get("blocks", [{}])[0]
        res = rep["parsed"].get("res") or {}
        cc = rep["checkcif"].get("counts") or {}
        cls = rep.get("checkcif_a_classification") or {}
        rc = rep.get("recompute") or {}
        rcc = rep.get("recompute_from_cif") or {}
        nomask = (rc.get("no_mask") or {}).get("R1_gt")
        rows.append({
            "rid": rid, "cond": man["condition"], "case": man["case_id"], "cond_zh": COND_ZH.get(man["condition"], man["condition"]),
            "r1_rep": v.get("r1_reported"), "r1_rec": v.get("r1_recomputed"), "r1_diff": v.get("r1_abs_difference"),
            "rec_route": "CIF 模型" if rcc else ("RES 零周期" if rc.get("main") else "—"),
            "nomask": nomask, "wr2": cif.get("wR2"), "goof": cif.get("GooF"), "dmin": cif.get("d_min"), "rint": cif.get("R_int"),
            "compl": cif.get("completeness"), "nrefl": cif.get("n_reflns_ls"), "flack": cif.get("flack"),
            "A": cc.get("A"), "A_subst": len(cls.get("substantive", [])), "B": cc.get("B"), "C": cc.get("C"),
            "status": v.get("scientific_status"), "lo": v.get("score_lower_bound"), "hi": v.get("score_upper_bound"), "unknown": v.get("unknown_items"),
            "tools": us.get("tool_calls"), "wall_min": round(man["science_seconds"] / 60, 1), "exit": man.get("exit_reason"),
            "early_stop": (man.get("early_stop") or {}).get("reason"), "mask": (res.get("flags") or {}).get("mask_abin"),
            "in_tok": us.get("input_tokens"), "out_tok": us.get("output_tokens"), "cost": us.get("cost_usd"),
            "claims_completion": (rep.get("claims") or {}).get("claims_completion"), "mentions_unresolved": (rep.get("claims") or {}).get("mentions_unresolved"),
        })
    L = []
    L.append("表 1　逐轮自动评价结果（12 轮）。R1 独立复算：带 SHELX 文件的轮次为 RES 零周期，纯产品轮次为 CIF 模型转 SHELXL、固定原子只精修标度；分数区间为预注册评分表的已确认下界至上界，其中两项对纯产品条件结构性不利。\n")
    L.append("| 轮 | 条件 | 案例 | R1 报告 | R1 独立复算 | 差值 | 复算途径 | 去掩膜 R1 | wR2 | 分辩率 Å | Rint | 完整度 | 精修反射 | checkCIF A（实质）/B/C | 科学状态 | 分数区间 | 工具调用 | 用时 min |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        L.append(f"| {r['rid']} | {r['cond_zh']} | {r['case']} | {fmt(r['r1_rep'])} | {fmt(r['r1_rec'])} | {fmt(r['r1_diff'])} | {r['rec_route']} | {fmt(r['nomask'])} | "
                 f"{fmt(r['wr2'], 3)} | {fmt(r['dmin'], 2)} | {fmt(r['rint'], 3)} | {pct(r['compl'])} | {fmt(r['nrefl'], 0)} | {r['A']}（{r['A_subst']}）/{r['B']}/{r['C']} | "
                 f"{STATUS_ZH.get(r['status'], r['status'])} | {fmt(r['lo'], 1)} 至 {fmt(r['hi'], 1)} | {r['tools']} | {r['wall_min']} |")
    L.append("\n表 2　按案例汇总\n")
    for case, title in (("nu1000", "NU-1000（给定 HKL，Zr-MOF，弱高角数据）"), ("alanine", "丙氨酸（646 帧 Rigaku 原始数据，从头还原）")):
        L.append(f"{title}\n")
        L.append("| 条件 | R1 复算 | 完整度 | 分辩率 | 实质 A 警报 | 科学状态 | 分数区间 | 用时 min | 工具 |")
        L.append("|---|---|---|---|---|---|---|---|---|")
        for r in rows:
            if r["case"] == case:
                L.append(f"| {r['cond_zh']} | {fmt(r['r1_rec'])} | {pct(r['compl'])} | {fmt(r['dmin'], 2)} | {r['A_subst']} | {STATUS_ZH.get(r['status'], r['status'])} | {fmt(r['lo'], 1)} 至 {fmt(r['hi'], 1)} | {r['wall_min']} | {r['tools']} |")
        L.append("")
    L.append("表 3　按条件汇总（两案例平均）\n")
    L.append("| 条件 | 平均分数区间 | 平均用时 min | 平均工具调用 | 实质 A 警报合计 | 输入 token 合计 | 输出 token 合计 |")
    L.append("|---|---|---|---|---|---|---|")
    for cond in ORDER:
        rs = [r for r in rows if r["cond"] == cond]
        lo = sum(r["lo"] for r in rs) / len(rs); hi = sum(r["hi"] for r in rs) / len(rs)
        L.append(f"| {COND_ZH[cond]} | {lo:.1f} 至 {hi:.1f} | {sum(r['wall_min'] for r in rs)/len(rs):.1f} | {sum(r['tools'] for r in rs)/len(rs):.0f} | {sum(r['A_subst'] for r in rs)} | "
                 f"{sum(r['in_tok'] or 0 for r in rs):,} | {sum(r['out_tok'] or 0 for r in rs):,} |")
    L.append("\n说明：input token 对 Claude Code 组不可比（其结果事件只报告未缓存的输入 token，几百个；实际花费见 usage.json 的 cost_usd，每轮 9.2 至 11.5 美元）。")
    L.append("\n表 4　自我判断与声明\n")
    L.append("| 轮 | 条件 | 案例 | 最终说明是否宣称完成 | 是否列出未解决问题 | 提前停止 |")
    L.append("|---|---|---|---|---|---|")
    for r in rows:
        L.append(f"| {r['rid']} | {r['cond_zh']} | {r['case']} | {'是' if r['claims_completion'] else '否'} | {'是' if r['mentions_unresolved'] else '否'} | {r['early_stop'] or '无'} |")
    (CTRL / "analysis" / "COMPARISON_12_tables.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    (CTRL / "analysis" / "comparison_rows.json").write_text(json.dumps(rows, ensure_ascii=False, indent=1), encoding="utf-8")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
