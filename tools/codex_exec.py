"""Supervised non-interactive run of the vendored Codex kernel (condition C0, and its smoke tests).

Launches `codex.exe exec --json ...` with the experiment's isolated CODEX_HOME, feeds the prompt on
stdin, streams the JSONL event log to disk line by line, enforces a hard wall-clock limit on the whole
process tree, and writes a commands/ record (argv, cwd, exit code, timing). Never uses `resume --last`.

Usage: codex_exec.py --cwd DIR --prompt FILE --out DIR --limit-seconds N [--image FILE ...]
"""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import psutil

ROOT = Path(r"H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921")
CODEX = ROOT / "engine" / "CrystalPilot" / "vendor" / "codex" / "node_modules" / "@openai" / "codex-win32-x64" \
    / "vendor" / "x86_64-pc-windows-msvc" / "bin" / "codex.exe"
CODEX_HOME = ROOT / "engine" / "codex-home-c0"
DIALS = Path(r"C:\users\<user>\miniforge3\envs\dials")
SHELX = ROOT / "toolbox" / "shelx"


def kill_tree(pid: int, events) -> None:
    try:
        root = psutil.Process(pid)
    except psutil.NoSuchProcess:
        return
    procs = root.children(recursive=True) + [root]
    for p in procs:
        try:
            events.write(json.dumps({"kind": "supervisor_kill", "pid": p.pid, "name": p.name(),
                                     "cmd": " ".join(p.cmdline())[:200], "ts": time.time()}) + "\n")
            p.kill()
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    psutil.wait_procs(procs, timeout=15)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cwd", required=True)
    ap.add_argument("--prompt", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--limit-seconds", type=int, required=True)
    ap.add_argument("--image", action="append", default=[])
    ap.add_argument("--model", default="gpt-6-astra")
    ap.add_argument("--effort", default="xhigh")
    ap.add_argument("--codex-home", default=str(CODEX_HOME))
    ap.add_argument("--sandbox", default="danger-full-access")
    a = ap.parse_args()

    out = Path(a.out)
    (out / "raw_events").mkdir(parents=True, exist_ok=True)
    (out / "commands").mkdir(parents=True, exist_ok=True)
    prompt = Path(a.prompt).read_text(encoding="utf-8")
    argv = [str(CODEX), "exec", "--json", "--sandbox", a.sandbox, "--skip-git-repo-check",
            "--model", a.model, "-c", f'model_reasoning_effort="{a.effort}"',
            "--output-last-message", str(out / "final_message.txt")]
    for img in a.image:
        argv += ["--image", img]
    argv.append("-")
    # a minimal, explicit environment: the kernel home, the science tools' directories, system dirs.
    env = {
        "CODEX_HOME": a.codex_home, "RUST_LOG": "warn", "PYTHONUTF8": "1",
        "SYSTEMROOT": os.environ.get("SYSTEMROOT", r"C:\Windows"), "COMSPEC": os.environ.get("COMSPEC", r"C:\Windows\System32\cmd.exe"),
        "TEMP": os.environ.get("TEMP", ""), "TMP": os.environ.get("TMP", ""),
        "USERPROFILE": os.environ.get("USERPROFILE", ""), "APPDATA": os.environ.get("APPDATA", ""),
        "LOCALAPPDATA": os.environ.get("LOCALAPPDATA", ""), "PROGRAMDATA": os.environ.get("PROGRAMDATA", ""),
        "HOMEDRIVE": os.environ.get("HOMEDRIVE", ""), "HOMEPATH": os.environ.get("HOMEPATH", ""),
        "PATHEXT": os.environ.get("PATHEXT", ".COM;.EXE;.BAT;.CMD;.PS1"),
        "PATH": ";".join([str(SHELX), str(DIALS / "Library" / "bin"), str(DIALS / "Scripts"), str(DIALS),
                          r"C:\Windows\System32", r"C:\Windows", r"C:\Windows\System32\WindowsPowerShell\v1.0",
                          r"C:\Windows\System32\Wbem"]),
        "NUMBER_OF_PROCESSORS": "4", "OMP_NUM_THREADS": "4",
    }
    rec = {"argv": [x for x in argv], "cwd": a.cwd, "env_keys": sorted(env), "path": env["PATH"],
           "started_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "limit_seconds": a.limit_seconds}
    (out / "commands" / "codex_exec.json").write_text(json.dumps(rec, ensure_ascii=False, indent=2), encoding="utf-8")

    t0 = time.monotonic()
    with (out / "raw_events" / "codex_events.jsonl").open("a", encoding="utf-8") as ev, \
            (out / "commands" / "codex_stderr.txt").open("a", encoding="utf-8") as err:
        proc = subprocess.Popen(argv, cwd=a.cwd, env=env, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                stderr=err, text=True, encoding="utf-8", errors="replace", bufsize=1)
        ev.write(json.dumps({"kind": "supervisor_start", "pid": proc.pid, "ts": time.time()}) + "\n"); ev.flush()
        try:
            proc.stdin.write(prompt)
            proc.stdin.close()
        except OSError:
            pass
        try:
            psutil.Process(proc.pid).nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
        except Exception:  # noqa: BLE001
            pass
        exit_reason = "completed"
        import threading
        lines = []

        def pump():
            for line in proc.stdout:
                ev.write(line if line.endswith("\n") else line + "\n"); ev.flush()
                lines.append(1)
        th = threading.Thread(target=pump, daemon=True)
        th.start()
        while True:
            if proc.poll() is not None:
                break
            if time.monotonic() - t0 > a.limit_seconds:
                exit_reason = "timeout"
                ev.write(json.dumps({"kind": "supervisor_timeout", "after_s": round(time.monotonic() - t0), "ts": time.time()}) + "\n"); ev.flush()
                kill_tree(proc.pid, ev)
                break
            time.sleep(2)
        th.join(timeout=30)
        # the sandbox helper or a science program may outlive the kernel: sweep the tree once more
        kill_tree(proc.pid, ev)
        rc = proc.poll()
    rec.update({"ended_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "wall_seconds": round(time.monotonic() - t0, 1),
                "exit_code": rc, "exit_reason": exit_reason, "event_lines": len(lines)})
    (out / "commands" / "codex_exec.json").write_text(json.dumps(rec, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({k: rec[k] for k in ("wall_seconds", "exit_code", "exit_reason", "event_lines")}))
    return 0 if exit_reason == "completed" else 3


if __name__ == "__main__":
    sys.exit(main())
