"""Aggregate sealed runs and verifier reports into control/analysis/: results.csv (one row per planned run,
including rows not run), run_index.csv (evidence paths), and summary.json for the report and demo page."""
import csv
import json
import time
from pathlib import Path

ROOT = Path(r"H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921")
CTRL = ROOT / "control"
OUT = CTRL / "analysis"
COLUMNS = ["run_id", "case_id", "condition", "replicate", "model_id", "kernel_version", "knowledge_mode", "run_status", "scientific_status",
           "publication_readiness", "protocol_validity", "start_utc", "end_utc", "wall_seconds", "agent_active_seconds", "exit_reason", "r1_reported", "r1_recomputed",
           "r1_definition", "wr2", "goof", "n_reflections", "resolution_angstrom", "reflection_set_hash", "mask_used", "excluded_reflections",
           "r1_abs_difference", "raw_pipeline_replay", "score_lower_bound", "score_upper_bound", "tool_calls", "model_requests", "input_tokens",
           "output_tokens", "cost_usd", "human_scientific_interventions", "automatic_interventions", "evidence_path", "plan_group", "concurrent_with", "early_stop_reason"]


def agent_active_seconds(evid: Path, cond: str) -> float | None:
    """Time from the task message to the agent's last turn end, from the native event stream. Differs from the
    runner's science_seconds when the runner waited after the agent had finished (after one rerun)."""
    p = evid / "raw_events" / "workbench_events.jsonl"
    if cond not in ("H0", "H1") or not p.exists():
        return None
    t_task, t_end = None, None
    for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
        try:
            rec = json.loads(line); e = rec["event"]
        except Exception:  # noqa: BLE001
            continue
        k = e.get("kind")
        if k == "user_message" and "inputs" in (e.get("text") or ""):
            t_task = rec.get("ts") if t_task is None else t_task
        if k in ("turn_completed", "turn_failed") and t_task is not None:
            t_end = rec.get("ts")
    return round(t_end - t_task, 1) if t_task and t_end else None


def to_utc(local: str | None) -> str | None:
    if not local:
        return None
    try:
        t = time.strptime(local[:19], "%Y-%m-%dT%H:%M:%S")
        return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(time.mktime(t)))
    except ValueError:
        return local


PLANS = [("plan.json", "pre-registered (frozen 2026-09-21 19:22)"),
         ("plan_extra_userrequest.json", "post-hoc rerun requested by the user after seeing r02/r03"),
         ("plan_bare.json", "bare product baselines 2026-09-21 22:20 (superseded 2026-09-22 01:40; r09 aborted, r10-r12 never run)"),
         ("plan_pure.json", "pure product baselines 2026-09-22 01:40 (no crystallographic software or library at all; plain Python only)")]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    plan = json.loads((CTRL / "freeze" / "plan.json").read_text(encoding="utf-8"))
    all_runs = []
    for fname, group in PLANS:
        p = CTRL / "freeze" / fname
        if p.exists():
            pl = json.loads(p.read_text(encoding="utf-8"))
            for r in pl["runs"]:
                all_runs.append({**r, "plan_group": group})
    state = json.loads((CTRL / "controller" / "controller_state.json").read_text(encoding="utf-8"))
    rows, index, summary = [], [], {"generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "runs": {}}
    for run in all_runs:
        rid = run["run_id"]
        st = state["runs"].get(rid, {})
        evid = CTRL / "runs" / rid
        man = json.loads((evid / "manifest.json").read_text(encoding="utf-8")) if (evid / "manifest.json").exists() else {}
        usage = json.loads((evid / "usage.json").read_text(encoding="utf-8")) if (evid / "usage.json").exists() else {}
        rep_p = CTRL / "verifier" / rid / "report.json"
        rep = json.loads(rep_p.read_text(encoding="utf-8")) if rep_p.exists() else {}
        v = rep.get("verdict", {})
        cif = ((rep.get("parsed") or {}).get("cif") or {}).get("blocks", [{}])
        cif = cif[0] if cif else {}
        res = (rep.get("parsed") or {}).get("res") or {}
        status = st.get("status") or run.get("status") or "planned"
        run_status = {"sealed": "completed", "sealed_protocol_invalid": "completed", "not_run_time": "not_run", "infra_failed": "infra_failed",
                      "planned": "planned", "running": "running", "preparing": "running", "collecting": "running",
                      "sealed_aborted": "aborted_by_user", "not_run_superseded": "not_run"}.get(status, status)
        if status in ("sealed", "sealed_protocol_invalid") and rep:
            run_status = "verified"
        exit_reason = man.get("exit_reason") or st.get("exit_reason") or st.get("reason")
        n_refl = cif.get("n_reflns_ls") or ((rep.get("recompute") or {}).get("main") or {}).get("n_all")
        row = {
            "run_id": rid, "case_id": run["case_id"], "condition": run["condition"], "replicate": run["replicate"],
            "model_id": "gpt-6-astra", "kernel_version": man.get("kernel_version") or plan["conditions"].get(run["condition"], {}).get("kernel"),
            "plan_group": run.get("plan_group"), "concurrent_with": man.get("concurrent_with"),
            "knowledge_mode": man.get("knowledge_mode"), "run_status": run_status,
            "scientific_status": v.get("scientific_status") if rep else ("not_evaluable" if run_status in ("not_run", "infra_failed") else None),
            "publication_readiness": v.get("publication_readiness") if rep else None,
            "protocol_validity": man.get("protocol_validity"),
            "start_utc": to_utc(man.get("started_science_at")), "end_utc": to_utc(man.get("ended_at")),
            "wall_seconds": man.get("science_seconds"), "agent_active_seconds": agent_active_seconds(evid, run["condition"]) if man else None, "exit_reason": exit_reason,
            "r1_reported": v.get("r1_reported"),
            "r1_recomputed": v.get("r1_recomputed") if v.get("r1_recomputed") is not None else (((rep.get("recompute_from_cif") or {}).get("closest_variant") or {}).get("R1_gt")),
            "r1_definition": v.get("r1_definition") if v.get("r1_recomputed") is not None else ("delivered CIF model + delivered reflections via SHELXL L.S. 0 (" + str(((rep.get("recompute_from_cif") or {}).get("closest_variant") or {}).get("variant")) + ")" if rep.get("recompute_from_cif") else v.get("r1_definition")),
            "wr2": cif.get("wR2"), "goof": cif.get("GooF"), "n_reflections": n_refl, "resolution_angstrom": cif.get("d_min"),
            "reflection_set_hash": ((rep.get("parsed") or {}).get("hkl") or {}).get("sha256"),
            "mask_used": (res.get("flags") or {}).get("mask_abin") if res else None,
            "excluded_reflections": (res.get("flags") or {}).get("omit") if res else None,
            "r1_abs_difference": v.get("r1_abs_difference"),
            "raw_pipeline_replay": (rep.get("chain") or {}).get("note") if run["case_id"] == "alanine" and rep else ("n/a (HKL provided)" if run["case_id"] == "nu1000" and rep else None),
            "score_lower_bound": v.get("score_lower_bound"), "score_upper_bound": v.get("score_upper_bound"),
            "tool_calls": usage.get("tool_calls"), "model_requests": usage.get("model_requests"),
            "input_tokens": usage.get("input_tokens"), "output_tokens": usage.get("output_tokens"), "cost_usd": usage.get("cost_usd"),
            "human_scientific_interventions": 0 if man else None, "automatic_interventions": man.get("interventions_auto"),
            "evidence_path": str(evid) if man else None,
            "early_stop_reason": (man.get("early_stop") or {}).get("reason") if man else None,
        }
        if run_status == "not_run":
            row["exit_reason"] = st.get("reason")
        rows.append(row)
        index.append({"run_id": rid, "arena": man.get("arena"), "evidence": str(evid) if man else None, "verifier": str(rep_p.parent) if rep else None,
                      "raw_events": str(next(iter((evid / "raw_events").glob("*.jsonl")), "")) if (evid / "raw_events").exists() else None,
                      "final_message": str(evid / "final_message.txt") if (evid / "final_message.txt").exists() else None,
                      "thread_or_session": man.get("thread_or_session")})
        summary["runs"][rid] = {"row": row, "verdict": v, "fatal": v.get("fatal_findings"), "claims": rep.get("claims"), "checkcif_counts": (rep.get("checkcif") or {}).get("counts"),
                                "delivery": {k: (rep.get("delivery") or {}).get(k) for k in ("formal_delivery", "picked")}, "final_text_tail": (evid / "final_message.txt").read_text(encoding="utf-8", errors="replace")[-1500:] if (evid / "final_message.txt").exists() else None}
    with (OUT / "results.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in COLUMNS})
    with (OUT / "run_index.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(index[0].keys()))
        w.writeheader()
        w.writerows(index)
    (OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=1, default=str), encoding="utf-8")
    for r in rows:
        print(f"{r['run_id']:18s} {r['run_status']:12s} {str(r['scientific_status']):24s} R1={r['r1_reported']} recomputed={r['r1_recomputed']} score={r['score_lower_bound']}-{r['score_upper_bound']} tools={r['tool_calls']} wall={r['wall_seconds']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
