"""Progress snapshot of the running pure-baseline trial(s), for the controller's no-progress checkpoints (early-stop-rule).

Prints: elapsed minutes, command / tool-use counts and their recent timing (idle detection), the agent's last
messages, what the agent has written into the arena (files by extension, newest first), and the rule audit
(crystallographic programs or libraries, other interpreters, downloads or installs). Read-only."""
import collections
import json
import re
import sys
import time
from pathlib import Path

ROOT = Path(r"H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921")
CTRL = ROOT / "control"
VIOL = re.compile(r"shelx|olex2|platon|dials|\bxds\b|crysalis|cctbx|gemmi|pymatgen|\base\b|spglib|crystalpilot(?!-campaigns)|miniforge|pyenv|conda|"
                  r"pip install|winget|npm install|Invoke-WebRequest|curl |wget |git clone|Start-BitsTransfer|D:\\OLEX2|H:\\CrystalPilot\\", re.I)
SKIP_DIRS = {"inputs", ".claude", "__pycache__"}


def running_runs() -> list:
    st = json.loads((CTRL / "controller" / "controller_state.json").read_text(encoding="utf-8"))
    return [rid for rid, s in st["runs"].items() if s.get("status") in ("running", "preparing", "collecting") and not rid.startswith("smoke")]


def commands_and_messages(evid: Path) -> tuple:
    cmds, msgs = [], []
    cs, cx = evid / "raw_events" / "claude_stream.jsonl", evid / "raw_events" / "codex_events.jsonl"
    if cs.exists():
        for l in cs.read_text(encoding="utf-8", errors="replace").splitlines():
            try:
                e = json.loads(l)
            except ValueError:
                continue
            if e.get("type") == "assistant":
                for b in (e.get("message") or {}).get("content") or []:
                    if not isinstance(b, dict):
                        continue
                    if b.get("type") == "tool_use":
                        inp = b.get("input") or {}
                        cmds.append((b.get("name"), str(inp.get("command") or inp.get("file_path") or json.dumps(inp, ensure_ascii=False))[:300]))
                    elif b.get("type") == "text" and b.get("text"):
                        msgs.append(b["text"])
    if cx.exists():
        for l in cx.read_text(encoding="utf-8", errors="replace").splitlines():
            try:
                e = json.loads(l)
            except ValueError:
                continue
            it = e.get("item") or {}
            if e.get("type") == "item.completed":
                if it.get("type") == "command_execution":
                    cmds.append(("shell", (it.get("command") or "")[:300] + f"  -> exit {it.get('exit_code')}"))
                elif it.get("type") == "agent_message":
                    msgs.append(it.get("text") or "")
    return cmds, msgs


def main() -> int:
    rids = sys.argv[1:] or running_runs()
    if not rids:
        print("no running run"); return 0
    for rid in rids:
        evid, arena = CTRL / "runs" / rid, ROOT / "arena" / rid
        man = json.loads((evid / "manifest.json").read_text(encoding="utf-8")) if (evid / "manifest.json").exists() else {}
        st = json.loads((CTRL / "controller" / "controller_state.json").read_text(encoding="utf-8"))["runs"].get(rid, {})
        started = st.get("started_at") or man.get("prepared_at")
        try:
            t0 = time.mktime(time.strptime(started[:19], "%Y-%m-%dT%H:%M:%S"))
            elapsed = (time.time() - t0) / 60
        except Exception:  # noqa: BLE001
            elapsed = float("nan")
        cmds, msgs = commands_and_messages(evid)
        raw = next(iter((evid / "raw_events").glob("*.jsonl")), None)
        last_event_age = (time.time() - raw.stat().st_mtime) / 60 if raw else float("nan")
        print(f"=== {rid} [{man.get('condition')}/{man.get('case_id')}] elapsed {elapsed:.1f} min of {man.get('limit_seconds', 0)//60} | "
              f"commands/tool uses {len(cmds)} | last event {last_event_age:.1f} min ago | exit {man.get('exit_reason')}")
        for m in msgs[-3:]:
            print("  MSG:", m[:400].replace("\n", " "))
        for name, c in cmds[-6:]:
            print(f"  {name}: {c[:180]}".replace("\n", " "))
        # what the agent produced
        files = []
        for p in arena.rglob("*"):
            if p.is_file() and not any(part in SKIP_DIRS for part in p.relative_to(arena).parts):
                files.append(p)
        by_ext = collections.Counter(p.suffix.lower() or "(none)" for p in files)
        newest = sorted(files, key=lambda p: p.stat().st_mtime, reverse=True)[:8]
        print(f"  files written: {len(files)} | by ext: {dict(by_ext.most_common(10))}")
        for p in newest:
            print(f"    {time.strftime('%H:%M:%S', time.localtime(p.stat().st_mtime))} {p.stat().st_size:>9d}  {p.relative_to(arena)}")
        # rule audit
        hits = [(n, c) for n, c in cmds if VIOL.search(c)]
        print(f"  RULE AUDIT: {len(hits)} suspicious command(s)")
        for n, c in hits[:6]:
            print(f"    !! {n}: {c[:200]}".replace("\n", " "))
        deliv = arena / "deliverables"
        if deliv.exists():
            print("  deliverables:", [p.name for p in deliv.iterdir()][:12])
    return 0


if __name__ == "__main__":
    sys.exit(main())
