"""Ask the runner to stop a native-kernel run early (pre-registered no-progress rule, early-stop-rule).

    stop_run.py <run_id> "<reason in one sentence>" [--checkpoint-minutes N]

Writes control/controller/stop_<run_id>.json; run_trial.py polls for it every 2 s, kills the kernel tree, records
exit_reason=stopped_by_controller with this marker in the manifest and seals the evidence as usual. The driver then
verifies the run and moves on. Nothing is killed here by name or by pid."""
import json
import sys
import time
from pathlib import Path

ROOT = Path(r"H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921")
CTRL = ROOT / "control"


def main() -> int:
    args = [a for a in sys.argv[1:]]
    cp = None
    if "--checkpoint-minutes" in args:
        i = args.index("--checkpoint-minutes"); cp = float(args[i + 1]); del args[i:i + 2]
    if len(args) < 2:
        print(__doc__); return 2
    run_id, reason = args[0], " ".join(args[1:])
    man_p = CTRL / "runs" / run_id / "manifest.json"
    if not man_p.exists():
        print("no such run evidence:", man_p); return 2
    man = json.loads(man_p.read_text(encoding="utf-8"))
    if man.get("exit_reason"):
        print(f"{run_id} already ended ({man['exit_reason']}); nothing to stop"); return 1
    marker = {"run_id": run_id, "decided_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "reason": reason,
              "rule": "early-stop-rule pre-registered no-progress rule (checkpoint at 30 min, then every 15 min)", "checkpoint_minutes": cp,
              "decided_by": "controller (Claude Code session acting as experiment controller), not the solving agent"}
    (CTRL / "controller" / f"stop_{run_id}.json").write_text(json.dumps(marker, ensure_ascii=False, indent=2), encoding="utf-8")
    with (CTRL / "controller" / "controller_events.jsonl").open("a", encoding="utf-8") as fh:
        fh.write(json.dumps({"ts": marker["decided_at"], "run_id": run_id, "event": "early_stop_requested", **marker}, ensure_ascii=False) + "\n")
    print("stop requested for", run_id, "->", reason)
    return 0


if __name__ == "__main__":
    sys.exit(main())
