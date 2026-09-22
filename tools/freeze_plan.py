"""Freeze the plan before any formal result exists: run list and order (seed 20260921), budgets, prompts,
condition settings, knowledge snapshot (copy + manifest), staging manifests, controller code hashes.

    freeze_plan.py [--include-p1 yes|no --p1-reason TEXT] [--smoke]
"""
import argparse
import hashlib
import json
import random
import shutil
import time
from pathlib import Path

ROOT = Path(r"H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921")
CTRL = ROOT / "control"
FREEZE = CTRL / "freeze"
KIT = FREEZE / "kit"
KNOWLEDGE = ROOT / "engine" / "CrystalPilot" / "knowledge"


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def manifest(root: Path) -> dict:
    return {p.relative_to(root).as_posix(): {"sha256": sha(p), "size": p.stat().st_size}
            for p in sorted(root.rglob("*")) if p.is_file()}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--include-p1", choices=["yes", "no"], required=True)
    ap.add_argument("--p1-reason", default="")
    a = ap.parse_args()
    if (FREEZE / "plan.json").exists():
        raise SystemExit("plan.json already frozen; refusing to overwrite")

    # prompts, verbatim from the kit
    (FREEZE / "prompts").mkdir(parents=True, exist_ok=True)
    for name in ("alanine.txt", "nu1000.txt", "common_runtime.md"):
        shutil.copy2(KIT / "prompts" / name, FREEZE / "prompts" / name)
    # knowledge snapshot (already redacted) frozen as a copy + manifest
    if (FREEZE / "knowledge_snapshot").exists():
        shutil.rmtree(FREEZE / "knowledge_snapshot")
    shutil.copytree(KNOWLEDGE, FREEZE / "knowledge_snapshot")
    km = manifest(KNOWLEDGE)
    (FREEZE / "knowledge_snapshot_manifest.json").write_text(json.dumps(km, indent=1), encoding="utf-8")

    cases = {
        "alanine": {"limit_minutes": 120, "prompt_file": "alanine.txt", "inputs": "inputs/alanine/ (raw Rigaku frames + acquisition metadata)"},
        "nu1000": {"limit_minutes": 180, "prompt_file": "nu1000.txt", "inputs": "inputs/start.hkl, inputs/start.ins, inputs/synthesis.png",
                   "scope_sentence": "只解析出临近发表级的框架结构，孔道内部的溶剂分子及客体无需处理。"},
    }
    conditions = {
        "C0": {"name": "native Codex kernel + same science software", "kernel": "codex-cli 0.155.0 (vendored)", "model": "gpt-6-astra", "effort": "xhigh",
               "provider": "crystalpilot", "sandbox": "danger-full-access (the Windows sandbox rejects every command in non-interactive mode; the workbench's auto mode likewise runs escalated commands unsandboxed)",
               "approval_policy": "never", "mcp": None, "instructions": "kernel base instructions only; no AGENTS.md", "knowledge": None},
        "H0": {"name": "CrystalPilot, knowledge_mode=tools_only", "kernel": "codex-cli 0.155.0 (same binary)", "model": "gpt-6-astra", "effort": "xhigh", "provider": "crystalpilot",
               "permission_mode": "auto", "settings": {"knowledge_mode": "tools_only", "subagents": "off", "enable_specialists": False, "allow_iucr_upload": False,
                                                        "model_override": "gpt-6-astra", "effort_override": "xhigh", "model_provider_override": "crystalpilot"},
               "instructions": "AGENTS.md tools-only template", "knowledge": None},
        "H1": {"name": "CrystalPilot, knowledge_mode=full, frozen redacted knowledge", "kernel": "codex-cli 0.155.0 (same binary)", "model": "gpt-6-astra", "effort": "xhigh", "provider": "crystalpilot",
               "permission_mode": "auto", "settings": {"knowledge_mode": "full", "subagents": "off", "enable_specialists": False, "allow_iucr_upload": False,
                                                        "model_override": "gpt-6-astra", "effort_override": "xhigh", "model_provider_override": "crystalpilot"},
               "instructions": "AGENTS.md full template", "knowledge": "frozen snapshot, 26 cards, 3 redacted (see provenance/knowledge_redactions.json)"},
        "P1": {"name": "Claude Code 2.1.261 + same science software (product comparison)", "kernel": "claude 2.1.261", "model": "gpt-6-astra via the gateway's Anthropic-compatible route",
               "effort": "gateway default; not controllable on that route (unverified)", "provider": "crystalpilot gateway /openai/v1/messages",
               "permission_mode": "dontAsk with Bash/Read/Write/Edit/Glob/Grep allowed; WebFetch/WebSearch/Task disallowed", "instructions": "product system prompt; no CLAUDE.md",
               "included": a.include_p1 == "yes", "exclusion_reason": a.p1_reason or None},
    }
    rng = random.Random(20260921)
    block_a = [("C0", "alanine"), ("C0", "nu1000"), ("H1", "alanine"), ("H1", "nu1000")]
    block_b = [("H0", "alanine"), ("H0", "nu1000")]
    block_c = [("P1", "alanine"), ("P1", "nu1000")]
    rng.shuffle(block_a); rng.shuffle(block_b); rng.shuffle(block_c)
    order = block_a + block_b + block_c
    runs = []
    for i, (cond, case) in enumerate(order, 1):
        rid = f"r{i:02d}_{cond}_{case}"
        status = "planned"
        if cond == "P1" and a.include_p1 != "yes":
            status = "preflight_blocked"
        runs.append({"run_id": rid, "order": i, "condition": cond, "case_id": case, "replicate": 1, "status": status,
                     "limit_minutes": cases[case]["limit_minutes"]})
    plan = {
        "frozen_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "seed": 20260921, "design": "serial, one run at a time, one replicate per condition x case",
        "why_serial": "one heavy science task at a time (DIALS full-frame integration, SHELXL on 337k reflections); machine has a history of memory-related crashes; ~24 h left",
        "drop_order_if_time_runs_out": "last rows first (P1, then H0); dropped rows stay in the results table as not run",
        "cases": cases, "conditions": conditions, "runs": runs,
        "budgets": {"alanine_minutes": 120, "nu1000_minutes": 180, "model_requests_soft": 200, "tool_calls_soft": 250,
                    "cost_note": "no reliable per-call price for the gateway; token counts are recorded; no automatic budget stop beyond wall clock"},
        "interventions_policy": {"approval_requests": "reject and log", "agent_question": "one fixed no-extra-information reply; for NU-1000 scope questions the verbatim scope sentence once",
                                 "background_job": "one fixed notice if the product does not surface it", "human_scientific": "none"},
        "hard_stop_local": "2026-09-22T20:00:00+08:00", "deadline_local": "2026-09-22T23:00:00+08:00",
    }
    (FREEZE / "plan.json").write_text(json.dumps(plan, ensure_ascii=False, indent=2), encoding="utf-8")
    frozen = {
        "frozen_at": plan["frozen_at"],
        "prompts": {n: sha(FREEZE / "prompts" / n) for n in ("alanine.txt", "nu1000.txt", "common_runtime.md")},
        "kit_files": manifest(KIT),
        "staging_manifests": {c: json.loads((ROOT / "staging" / c / "manifest.json").read_text(encoding="utf-8")) for c in ("alanine", "nu1000")},
        "knowledge_snapshot_files": len(km), "knowledge_snapshot_sha256_of_manifest": hashlib.sha256(json.dumps(km, sort_keys=True).encode()).hexdigest(),
        "controller_code": {p.name: sha(p) for p in sorted((CTRL / "tools").glob("*.py"))},
        "kernel_homes": {h: {f.name: sha(f) for f in sorted((ROOT / "engine" / h).glob("*")) if f.is_file()} for h in ("codex-home-wb", "codex-home-c0")},
        "engine_commit": "327e1ef7bcb3af7c4922dca1cf83ef7b1af63cdf",
        "grader": "03_独立评分与失败分析.md rubric; verifier scripts in control/tools/verify_run.py (written before the first result is scored)",
    }
    (FREEZE / "protocol_frozen.json").write_text(json.dumps(frozen, ensure_ascii=False, indent=2), encoding="utf-8")
    print("frozen:", plan["frozen_at"])
    for r in runs:
        print(f"  {r['run_id']:22s} {r['status']:18s} {r['limit_minutes']} min")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
