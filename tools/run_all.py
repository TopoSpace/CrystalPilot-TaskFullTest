"""Unattended driver: runs the frozen plan in order, one trial at a time, verifies each sealed run,
and stops at the first run that cannot finish before the science hard stop (later rows are marked not run,
which is the pre-registered drop order). Meant to be started detached; progress goes to driver.log and the
controller state file.
"""
import json
import subprocess
import sys
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

ROOT = Path(r"H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921")
CTRL = ROOT / "control"
ENGINE = ROOT / "engine" / "CrystalPilot"
PY = ENGINE / ".venv" / "Scripts" / "python.exe"
STATE = CTRL / "controller" / "controller_state.json"
LOG = CTRL / "controller" / "driver.log"
HARD_STOP = datetime(2026, 9, 22, 20, 0, tzinfo=timezone(timedelta(hours=8)))
MARGIN_S = 25 * 60


def log(msg: str) -> None:
    line = f"[{time.strftime('%Y-%m-%dT%H:%M:%S%z')}] {msg}"
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


def set_run(run_id: str, **fields) -> None:
    st = json.loads(STATE.read_text(encoding="utf-8"))
    st["runs"].setdefault(run_id, {}).update(fields)
    st["updated_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    STATE.write_text(json.dumps(st, ensure_ascii=False, indent=2), encoding="utf-8")


def server_ok() -> bool:
    import urllib.request
    try:
        with urllib.request.urlopen("http://127.0.0.1:8021/api/health", timeout=10) as r:
            return r.status == 200
    except Exception:  # noqa: BLE001
        return False


sys.path.insert(0, str(Path(__file__).resolve().parent))
import cpucap  # noqa: E402  (same directory)


def main() -> int:
    args = [a for a in sys.argv[1:]]
    plan_path = CTRL / "freeze" / "plan.json"
    if "--plan" in args:
        i = args.index("--plan"); plan_path = Path(args[i + 1]); del args[i:i + 2]
    log("cpu cap: " + cpucap.limit_self())
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    start_from = args[0] if args else None
    st = json.loads(STATE.read_text(encoding="utf-8"))
    st["phase"] = "run"
    st["driver_started_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    STATE.write_text(json.dumps(st, ensure_ascii=False, indent=2), encoding="utf-8")
    log(f"driver start; plan {plan_path.name} frozen {plan.get('amendment', {}).get('frozen_at') or plan['frozen_at']}; {len(plan['runs'])} runs")
    active = start_from is None
    for run in plan["runs"]:
        rid = run["run_id"]
        if not active:
            active = rid == start_from
            if not active:
                continue
        st = json.loads(STATE.read_text(encoding="utf-8"))
        if st["runs"].get(rid, {}).get("status") in ("sealed", "sealed_protocol_invalid", "verified", "not_run_time"):
            log(f"{rid}: already {st['runs'][rid]['status']}, skipping")
            continue
        now = datetime.now(timezone(timedelta(hours=8)))
        need = run["limit_minutes"] * 60 + MARGIN_S
        remaining = (HARD_STOP - now).total_seconds()
        if need > remaining:
            log(f"{rid}: needs {need/60:.0f} min but only {remaining/60:.0f} min remain before the hard stop; this and all later rows are not run")
            for later in plan["runs"][plan["runs"].index(run):]:
                set_run(later["run_id"], status="not_run_time", reason="insufficient time before the science hard stop (pre-registered drop order: last rows first)")
            break
        if run["condition"] in ("H0", "H1") and not server_ok():
            log(f"{rid}: workbench not healthy; starting it")
            subprocess.run(["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", str(CTRL / "tools" / "start_server.ps1")],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=300)
            if not server_ok():
                log(f"{rid}: workbench still down; marking infra_failed and continuing")
                set_run(rid, status="infra_failed", reason="workbench server could not be started")
                continue
        log(f"{rid}: starting ({run['condition']} / {run['case_id']} / {run['limit_minutes']} min)")
        t0 = time.time()
        with (CTRL / "controller" / f"driver_{rid}.log").open("a", encoding="utf-8") as out:
            rc = subprocess.run([str(PY), "-X", "utf8", str(CTRL / "tools" / "run_trial.py"), "--run-id", rid, "--plan", str(plan_path)],
                                cwd=str(ENGINE), stdout=out, stderr=subprocess.STDOUT, timeout=(run["limit_minutes"] + 60) * 60).returncode
        log(f"{rid}: run_trial exit {rc} after {(time.time() - t0)/60:.1f} min")
        # independent verification of the sealed run (read-only on arena; failures are logged, not fatal)
        t1 = time.time()
        try:
            with (CTRL / "controller" / f"driver_{rid}.log").open("a", encoding="utf-8") as out:
                vrc = subprocess.run([str(PY), "-X", "utf8", str(CTRL / "tools" / "verify_run.py"), "--run-id", rid],
                                     cwd=str(ENGINE), stdout=out, stderr=subprocess.STDOUT, timeout=2400).returncode
            log(f"{rid}: verify_run exit {vrc} after {(time.time() - t1)/60:.1f} min")
            if vrc == 0:
                set_run(rid, verified_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
        except subprocess.TimeoutExpired:
            log(f"{rid}: verify_run timed out")
        time.sleep(20)
    st = json.loads(STATE.read_text(encoding="utf-8"))
    st["phase"] = "runs_finished"
    st["driver_finished_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    STATE.write_text(json.dumps(st, ensure_ascii=False, indent=2), encoding="utf-8")
    log("driver finished")
    return 0


if __name__ == "__main__":
    sys.exit(main())
