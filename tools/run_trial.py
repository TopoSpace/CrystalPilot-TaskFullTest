"""Run one pre-registered trial unattended and seal its evidence.

    run_trial.py --run-id <id> [--limit-minutes N]   (the limit override exists for the mechanics smoke test only)

Reads control/freeze/plan.json, prepares arena/<run_id>/ from staging, starts the condition's backend
(C0 native Codex, H0/H1 CrystalPilot workbench, P1 Claude Code), supervises it against the hard wall-clock
limit, applies only the fixed-template interventions the protocol allows, then collects and seals.
Every state change is written to control/controller/controller_state.json and controller_events.jsonl.
"""
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

import psutil

sys.path.insert(0, str(Path(__file__).parent))
from cp_client import Sink, Workbench  # noqa: E402
import cpucap  # noqa: E402

ROOT = Path(r"H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921")
CTRL = ROOT / "control"
FREEZE = CTRL / "freeze"
STATE = CTRL / "controller" / "controller_state.json"
EVENTS = CTRL / "controller" / "controller_events.jsonl"
LOCK = CTRL / "controller" / "active_run.lock"
DIALS = Path(r"C:\users\<user>\miniforge3\envs\dials")
TOOLBOX = ROOT / "toolbox"
CODEX = ROOT / "engine/CrystalPilot/vendor/codex/node_modules/@openai/codex-win32-x64/vendor/x86_64-pc-windows-msvc/bin/codex.exe"
CODEX_HOME_C0 = ROOT / "engine" / "codex-home-c0"
CLAUDE_HOME = ROOT / "engine" / "claude-home"
KNOWLEDGE = ROOT / "engine" / "CrystalPilot" / "knowledge"
HARD_STOP = "2026-09-22T20:00:00+0800"

WARMUP = "不开始任何任务，也不要读取任何文件。只回复一行“就绪”。"
REPLY_NO_INFO = "没有额外信息，请基于已提供材料继续；无法确定的内容请标为未知并说明影响。"
REPLY_SCOPE = "只解析出临近发表级的框架结构，孔道内部的溶剂分子及客体无需处理。"
NOTICE_BG = "后台作业已结束/仍在运行，请在本轮预算内检查并完成交付。"
QUESTION_MARKERS = ("是否", "需要我", "请确认", "要不要", "可以吗", "请告知", "请提供", "should i", "do you want", "please confirm", "would you like", "let me know")
SCOPE_MARKERS = ("孔道", "客体", "溶剂", "guest", "solvent", "pore")

#: bare product baselines (C0B / P1B, added 2026-09-21 22:20 at the user's request): no software inventory, the
#: machine's ordinary PATH, and an explicit ban on downloading or installing anything
BARE_PREFIX = ("环境约束：不允许下载、安装或更新任何软件、Python 包、模型或数据，只能使用本机现有的程序与文件；"
               "本轮不提供软件清单，需要什么工具请自行在本机查找。")
BARE_RULES = [
    "不允许下载、安装或更新任何软件、Python 包、模型或数据；只能使用本机现有的程序与文件。",
    "本文件不列出软件清单；需要什么工具请自行在本机查找。",
    "不访问网络查找该样品的答案。",
]

#: pure product baselines (C0N / P1N, user redesign 2026-09-22 01:40): no crystallographic program or library of
#: any kind may be used; only a plain Python (numpy/scipy/matplotlib/pandas/pillow) and Windows built-ins
PURE_PY = TOOLBOX / "pyplain"
PURE_PREFIX = ("环境约束（本轮特别规定）：本轮不提供、也不允许使用任何晶体学软件或晶体学程序库，包括但不限于 SHELX 系列"
               "（SHELXT、SHELXL、SHELXS、SHELXD 等）、Olex2、PLATON、DIALS、XDS、CrysAlisPro、cctbx、gemmi、pymatgen、ASE、spglib、"
               "CrystalPilot，以及本机上任何其他现成的数据还原、结构求解、精修或验证程序；也不允许下载、安装或更新任何软件、Python 包、"
               "模型或数据。你只能使用 environment.json 指定的通用 Python 环境（Python 3.12，仅含 numpy、scipy、matplotlib、pandas、pillow）"
               "和 Windows 自带的命令行工具，自行编写全部算法完成任务：读取输入数据、必要时的指标化与积分、结构求解、精修与检验都由你自己实现。"
               "不得调用本机其他 Python 解释器或其已安装的包，不得读取或复用本机任何现成程序的源代码。")
PURE_RULES = [
    "不允许使用任何晶体学软件或晶体学程序库（SHELX 系列、Olex2、PLATON、DIALS、XDS、CrysAlisPro、cctbx、gemmi、pymatgen、ASE、spglib、CrystalPilot 及本机其他现成的还原、求解、精修、验证程序）。",
    "不允许下载、安装或更新任何软件、Python 包、模型或数据。",
    "只能使用下面 environment 中列出的 Python 环境与 Windows 自带命令行工具；全部算法自行编写。",
    "不得调用本机其他 Python 解释器或其已安装的包，不得读取或复用本机任何现成程序的源代码。",
    "不访问网络查找该样品的答案。",
]


def pure_path() -> str:
    return ";".join([str(PURE_PY / "Scripts"), r"C:\Windows\System32", r"C:\Windows",
                     r"C:\Windows\System32\WindowsPowerShell\v1.0", r"C:\Windows\System32\Wbem"])


def pure_env_extra() -> dict:
    return {"PYTHONNOUSERSITE": "1", "VIRTUAL_ENV": str(PURE_PY), "PYTHONIOENCODING": "utf-8", "PYTHONUTF8": "1"}


def machine_path() -> str:
    """The PATH a normal shell on this machine has (system + user registry values), with nothing added by the
    experiment. Entries pointing into the experiment, the product checkout or its tool environments are dropped."""
    import winreg
    parts: list = []
    for hive, key in ((winreg.HKEY_LOCAL_MACHINE, r"SYSTEM\CurrentControlSet\Control\Session Manager\Environment"),
                      (winreg.HKEY_CURRENT_USER, r"Environment")):
        try:
            with winreg.OpenKey(hive, key) as k:
                val, _ = winreg.QueryValueEx(k, "Path")
                parts += [os.path.expandvars(p) for p in str(val).split(";") if p.strip()]
        except OSError:
            continue
    seen, out = set(), []
    for p in parts:
        low = p.lower().rstrip("\\")
        if any(s in low for s in ("crystalpilot", "scxrd-agent-eval", "envs\\dials", "vendor\\shelx")):
            continue
        if low not in seen:
            seen.add(low)
            out.append(p)
    return ";".join(out)


CONDITION_SETTINGS = {
    "subagents": "off", "enable_specialists": False, "allow_iucr_upload": False,
    "model_override": "gpt-6-astra", "effort_override": "xhigh", "model_provider_override": "crystalpilot",
}


# ----------------------------------------------------------------------------- helpers
def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def sha256_of(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def log_event(run_id: str, event: str, **extra) -> None:
    rec = {"ts": now(), "run_id": run_id, "event": event, **extra}
    with EVENTS.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec, ensure_ascii=False, default=str) + "\n")
    print(f"[{rec['ts']}] {run_id}: {event} {json.dumps(extra, ensure_ascii=False, default=str)[:200]}", flush=True)


def update_state(run_id: str, **fields) -> None:
    st = json.loads(STATE.read_text(encoding="utf-8"))
    run = st["runs"].setdefault(run_id, {})
    run.update(fields)
    run["updated_at"] = now()
    if not run.get("concurrent_with") and not fields.get("concurrent_with"):
        st["active_run"] = run_id if fields.get("status") in ("preparing", "warming_up", "running", "collecting") else (None if st.get("active_run") == run_id else st.get("active_run"))
    st["updated_at"] = now()
    STATE.write_text(json.dumps(st, ensure_ascii=False, indent=2), encoding="utf-8")


STRONG_ASKS = ("请补充", "请提供", "请告知", "请确认", "请说明您", "请给出", "需要您", "等待您", "please provide", "please confirm", "please let me know")


def looks_like_question(text: str) -> bool:
    t = (text or "").strip().lower()
    tail = t[-600:]
    if any(m in tail for m in STRONG_ASKS):
        return True
    return ("?" in tail or "？" in tail) and any(m in tail for m in QUESTION_MARKERS)


def mentions_scope(text: str) -> bool:
    t = (text or "").lower()
    return any(m in t for m in SCOPE_MARKERS)


def kill_tree(pid: int, sink) -> None:
    try:
        root = psutil.Process(pid)
    except psutil.NoSuchProcess:
        return
    procs = root.children(recursive=True) + [root]
    for p in procs:
        try:
            sink(json.dumps({"kind": "supervisor_kill", "pid": p.pid, "name": p.name(), "cmd": " ".join(p.cmdline())[:200], "ts": time.time()}))
            p.kill()
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    psutil.wait_procs(procs, timeout=20)


def stop_marker(run_id: str) -> dict | None:
    """Controller early-stop request for this run (written by stop_run.py), or None."""
    p = CTRL / "controller" / f"stop_{run_id}.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {"reason": "unreadable marker"}


class Interventions:
    def __init__(self, path: Path):
        self.path = path
        self.count_auto = 0
        self.count_human_scientific = 0

    def add(self, kind: str, **extra) -> None:
        self.count_auto += 1
        with self.path.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps({"ts": now(), "kind": kind, "automatic": True, **extra}, ensure_ascii=False) + "\n")


def manifest_dir(root: Path, skip_dirs=()) -> dict:
    out = {}
    for folder, dirs, files in os.walk(root):
        dirs[:] = [d for d in dirs if d not in skip_dirs]
        for f in files:
            p = Path(folder) / f
            rel = p.relative_to(root).as_posix()
            try:
                out[rel] = {"sha256": sha256_of(p), "size": p.stat().st_size}
            except OSError as e:
                out[rel] = {"error": type(e).__name__}
    return out


# ----------------------------------------------------------------------------- preparation
def prepare(run: dict, plan: dict, limit_s: int) -> dict:
    run_id, case, cond = run["run_id"], run["case_id"], run["condition"]
    arena = ROOT / "arena" / run_id
    evid = CTRL / "runs" / run_id
    if arena.exists():
        raise SystemExit(f"arena dir already exists for {run_id}; a rerun is a new trial id")
    arena.mkdir(parents=True)
    for sub in ("raw_events", "commands", "effective_instructions", "artifacts", "native_state"):
        (evid / sub).mkdir(parents=True, exist_ok=True)
    # inputs from the sealed staging copy; hashes re-checked against the staging manifest
    src = ROOT / "staging" / case / "inputs"
    shutil.copytree(src, arena / "inputs")
    staged = json.loads((ROOT / "staging" / case / "manifest.json").read_text(encoding="utf-8"))
    mismatch = []
    for rel, meta in staged.items():          # keys are relative to inputs/
        p = arena / "inputs" / rel
        if not p.exists() or sha256_of(p) != meta["sha256"]:
            mismatch.append(rel)
    if mismatch:
        raise SystemExit(f"input copy mismatch: {mismatch[:5]}")
    for p in (arena / "inputs").rglob("*"):
        if p.is_file():
            try:
                os.chmod(p, 0o444)
            except OSError:
                pass
    (arena / "deliverables").mkdir(exist_ok=True)
    bare = cond in ("C0B", "P1B")
    pure = cond in ("C0N", "P1N")
    common = (FREEZE / "prompts" / "common_runtime.md").read_text(encoding="utf-8").strip()
    task = (FREEZE / "prompts" / f"{plan['cases'][case]['prompt_file']}").read_text(encoding="utf-8").strip()
    prefix = PURE_PREFIX if pure else (BARE_PREFIX if bare else "")
    prompt = (prefix + "\n\n" if prefix else "") + common + "\n\n" + task
    (arena / "TASK.md").write_text(prompt + "\n", encoding="utf-8")
    (evid / "prompt.txt").write_text(prompt + "\n", encoding="utf-8")
    if pure:
        env_card = {
            "run_id": run_id, "case": case, "condition": cond,
            "budget": {"wall_clock_minutes": limit_s // 60, "note": "the controller interrupts the run at the wall-clock limit; deliver a checkable state before then"},
            "working_directory": str(arena), "inputs": str(arena / "inputs"), "inputs_are_read_only": True,
            "deliverables_directory": str(arena / "deliverables"),
            "deliverables_note": "put the final structure files (CIF and any other formats you produce), your code, the computation logs and a short written summary here",
            "rules": PURE_RULES,
            "environment": {
                "python": str(PURE_PY / "Scripts" / "python.exe"),
                "packages": ["numpy", "scipy", "matplotlib", "pandas", "pillow"],
                "shell": "Windows built-in cmd / PowerShell",
                "note": "this is the only Python you may use; it has no crystallographic library; PATH contains only this Python and the Windows system directories",
            },
            "software": "no crystallographic program or library is provided or may be used; implement every algorithm yourself",
            "language": "reply in the language of the task text",
        }
    elif bare:
        env_card = {
            "run_id": run_id, "case": case, "condition": cond,
            "budget": {"wall_clock_minutes": limit_s // 60, "note": "the controller interrupts the run at the wall-clock limit; deliver a checkable state before then"},
            "working_directory": str(arena), "inputs": str(arena / "inputs"), "inputs_are_read_only": True,
            "deliverables_directory": str(arena / "deliverables"),
            "deliverables_note": "put the final structure files, the computation inputs/outputs you consider final and a short written summary here",
            "rules": BARE_RULES,
            "software": "not listed in this run; nothing may be downloaded or installed",
            "language": "reply in the language of the task text",
        }
    else:
        env_card = {
        "run_id": run_id, "case": case, "condition": cond,
        "budget": {"wall_clock_minutes": limit_s // 60, "note": "the controller interrupts the run at the wall-clock limit; deliver a checkable state before then",
                   "model_requests_soft_limit": 200, "tool_calls_soft_limit": 250},
        "working_directory": str(arena), "inputs": str(arena / "inputs"), "inputs_are_read_only": True,
        "deliverables_directory": str(arena / "deliverables"),
        "deliverables_note": "put the final structure files (CIF/RES/INS/HKL/FCF as applicable), the computation inputs/outputs you consider final and a short written summary here; CrystalPilot users may also rely on the product's own results folder",
        "software": {
            "shelxl": str(TOOLBOX / "shelx" / "shelxl.exe"), "shelxt": str(TOOLBOX / "shelx" / "shelxt.exe"),
            "platon": str(TOOLBOX / "shelx" / "platon.exe"), "olex2": str(TOOLBOX / "olex2" / "app" / "olex2.exe"),
            "dials": {"prefix": str(DIALS), "dispatchers": str(DIALS / "Scripts"), "dll_dir": str(DIALS / "Library" / "bin"),
                      "usage": "run the dials.* dispatchers from Scripts with Library\\bin on PATH (the C0 shell already has both on PATH); dials.import understands Rigaku .rodhypix frames",
                      "version": "DIALS 3.30"},
            "python": {"path": str(DIALS / "python.exe"), "packages": ["cctbx", "iotbx", "smtbx", "dxtbx", "dials", "gemmi", "numpy", "scipy", "matplotlib", "PIL"],
                       "note": "a plain scientific Python; import errors for compiled extensions disappear when Library\\bin is on PATH"},
            "help": "each program prints its usage with -h/--help or when started without arguments",
        },
        "network": "no internet lookup of this sample; documentation of the software may be consulted",
        "language": "reply in the language of the task text",
        }
    (arena / "environment.json").write_text(json.dumps(env_card, ensure_ascii=False, indent=2), encoding="utf-8")
    manifest = {
        "run_id": run_id, "case_id": case, "condition": cond, "replicate": run["replicate"],
        "model_id": "gpt-6-astra", "reasoning_effort": "xhigh" if cond not in ("P1", "P1B", "P1N") else "gateway default (unverified)",
        "provider_id": "crystalpilot", "kernel_version": "codex-cli 0.155.0" if cond not in ("P1", "P1B", "P1N") else "claude 2.1.261",
        "knowledge_mode": {"C0": None, "C0B": None, "C0N": None, "H0": "tools_only", "H1": "full", "P1": None, "P1B": None, "P1N": None}[cond],
        "bare_baseline": bare, "pure_baseline": pure, "message_variant": "pure (no-tools constraint prepended)" if pure else ("bare (constraint sentence prepended)" if bare else "standard"),
        "limit_seconds": limit_s, "prompt_sha256": hashlib.sha256(prompt.encode("utf-8")).hexdigest(),
        "inputs_manifest": staged, "arena": str(arena), "prepared_at": now(),
    }
    (evid / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return {"arena": arena, "evid": evid, "prompt": prompt, "manifest": manifest}


# ----------------------------------------------------------------------------- backends
def run_c0(ctx: dict, run: dict, limit_s: int, bare: bool = False, pure: bool = False) -> dict:
    arena, evid = ctx["arena"], ctx["evid"]
    early_stop: dict = {}
    iv = Interventions(evid / "interventions.jsonl")
    ev_path = evid / "raw_events" / "codex_events.jsonl"
    err_path = evid / "commands" / "codex_stderr.txt"
    images = [str(arena / "inputs" / "synthesis.png")] if run["case_id"] == "nu1000" else []
    tool_path = pure_path() if pure else machine_path() if bare else ";".join([
        str(TOOLBOX / "shelx"), str(DIALS / "Library" / "bin"), str(DIALS / "Scripts"), str(DIALS),
        r"C:\Windows\System32", r"C:\Windows", r"C:\Windows\System32\WindowsPowerShell\v1.0", r"C:\Windows\System32\Wbem"])
    env = {
        "CODEX_HOME": str(CODEX_HOME_C0), "RUST_LOG": "warn", "PYTHONUTF8": "1",
        "SYSTEMROOT": os.environ.get("SYSTEMROOT", r"C:\Windows"), "COMSPEC": os.environ.get("COMSPEC", r"C:\Windows\System32\cmd.exe"),
        "TEMP": os.environ.get("TEMP", ""), "TMP": os.environ.get("TMP", ""), "USERPROFILE": os.environ.get("USERPROFILE", ""),
        "APPDATA": os.environ.get("APPDATA", ""), "LOCALAPPDATA": os.environ.get("LOCALAPPDATA", ""), "PROGRAMDATA": os.environ.get("PROGRAMDATA", ""),
        "HOMEDRIVE": os.environ.get("HOMEDRIVE", ""), "HOMEPATH": os.environ.get("HOMEPATH", ""), "PATHEXT": os.environ.get("PATHEXT", ".COM;.EXE;.BAT;.CMD;.PS1"),
        "PATH": tool_path,
        "NUMBER_OF_PROCESSORS": "4", "OMP_NUM_THREADS": "4",
    }
    if pure:
        env.update(pure_env_extra())
    (evid / "commands" / "environment_path.txt").write_text(tool_path.replace(";", "\n") + "\n", encoding="utf-8")
    base = [str(CODEX), "exec", "--json", "--sandbox", "danger-full-access", "--skip-git-repo-check",
            "--model", "gpt-6-astra", "-c", 'model_reasoning_effort="xhigh"']
    t_start = time.monotonic()
    deadline = t_start + limit_s
    thread_id = None
    last_text = ""
    replies_sent = {"no_info": 0, "scope": 0}
    attempts = []

    def one_attempt(argv, stdin_text, tag):
        nonlocal thread_id, last_text
        rec = {"tag": tag, "argv": argv, "cwd": str(arena), "started_at": now()}
        with ev_path.open("a", encoding="utf-8") as ev, err_path.open("a", encoding="utf-8") as err:
            proc = subprocess.Popen(argv, cwd=str(arena), env=env, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=err,
                                    text=True, encoding="utf-8", errors="replace", bufsize=1)
            try:
                psutil.Process(proc.pid).nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            except Exception:  # noqa: BLE001
                pass
            ev.write(json.dumps({"kind": "supervisor_start", "tag": tag, "pid": proc.pid, "ts": time.time()}) + "\n"); ev.flush()
            try:
                proc.stdin.write(stdin_text); proc.stdin.close()
            except OSError:
                pass
            reason = "completed"

            def pump():
                nonlocal thread_id, last_text
                for line in proc.stdout:
                    ev.write(line if line.endswith("\n") else line + "\n"); ev.flush()
                    try:
                        e = json.loads(line)
                    except ValueError:
                        continue
                    if e.get("type") == "thread.started":
                        thread_id = e.get("thread_id")
                    it = e.get("item") or {}
                    if e.get("type") == "item.completed" and it.get("type") == "agent_message":
                        last_text = it.get("text") or ""
            th = threading.Thread(target=pump, daemon=True); th.start()
            while proc.poll() is None:
                if time.monotonic() > deadline:
                    reason = "timeout"
                    ev.write(json.dumps({"kind": "supervisor_timeout", "ts": time.time()}) + "\n"); ev.flush()
                    kill_tree(proc.pid, lambda s: ev.write(s + "\n"))
                    break
                mk = stop_marker(run["run_id"])
                if mk:
                    reason = "stopped_by_controller"
                    early_stop.update(mk)
                    ev.write(json.dumps({"kind": "supervisor_stop", "ts": time.time(), "marker": mk}, ensure_ascii=False) + "\n"); ev.flush()
                    kill_tree(proc.pid, lambda s: ev.write(s + "\n"))
                    break
                time.sleep(2)
            th.join(timeout=30)
            kill_tree(proc.pid, lambda s: ev.write(s + "\n"))
            rec.update({"ended_at": now(), "exit_code": proc.poll(), "exit_reason": reason})
        attempts.append(rec)
        return reason

    reason = one_attempt(base + ["--output-last-message", str(evid / "final_message.txt")] + sum([["--image", i] for i in images], []) + ["-"], ctx["prompt"], "initial")
    # the protocol's single fixed replies, only if the agent stopped with a question and budget remains
    for _ in range(2):
        if reason != "completed" or time.monotonic() > deadline - 600 or not thread_id:
            break
        if not looks_like_question(last_text):
            break
        if run["case_id"] == "nu1000" and mentions_scope(last_text) and replies_sent["scope"] == 0:
            reply, key = REPLY_SCOPE, "scope"
        elif replies_sent["no_info"] == 0:
            reply, key = REPLY_NO_INFO, "no_info"
        else:
            break
        replies_sent[key] += 1
        iv.add("fixed_reply", template=key, text=reply, trigger=last_text[-300:])
        reason = one_attempt(base + ["--output-last-message", str(evid / f"final_message_{key}.txt"), "resume", thread_id, reply], "", f"resume_{key}")
    (evid / "commands" / "codex_attempts.json").write_text(json.dumps(attempts, ensure_ascii=False, indent=2), encoding="utf-8")
    return {"exit_reason": reason, "thread_id": thread_id, "science_seconds": round(time.monotonic() - t_start, 1),
            "warmup_seconds": 0, "interventions_auto": iv.count_auto, "final_text": last_text[-2000:], "early_stop": early_stop or None}


def run_p1(ctx: dict, run: dict, limit_s: int, bare: bool = False, pure: bool = False) -> dict:
    arena, evid = ctx["arena"], ctx["evid"]
    early_stop: dict = {}
    iv = Interventions(evid / "interventions.jsonl")
    stream = evid / "raw_events" / "claude_stream.jsonl"
    err_path = evid / "commands" / "claude_stderr.txt"
    # a clean environment: no marker of the controller's own agent session (a nested Claude Code otherwise
    # treats provider routing as host-managed and ignores the isolated settings), no inherited API variables
    env = {k: v for k, v in os.environ.items()
           if not k.upper().startswith(("CLAUDE", "ANTHROPIC", "MCP_", "OTEL", "CODEX", "CRYSTALPILOT", "OPENAI"))}
    env["CLAUDE_CONFIG_DIR"] = str(CLAUDE_HOME)
    if pure:
        env = {k: v for k, v in env.items() if not k.upper().startswith(("PYTHON", "VIRTUAL_ENV", "CONDA"))}
        env["PATH"] = pure_path()
        env.update(pure_env_extra())
    elif bare:
        env["PATH"] = machine_path()
    else:
        env["PATH"] = ";".join([str(TOOLBOX / "shelx"), str(DIALS / "Library" / "bin"), str(DIALS / "Scripts"), str(DIALS), env.get("PATH", "")])
    (evid / "commands" / "environment_path.txt").write_text(env["PATH"].replace(";", "\n") + "\n", encoding="utf-8")
    claude_cmd = os.environ.get("CLAUDE_CMD_PATH", r"C:\users\<user>\AppData\Roaming\npm\claude.CMD")
    base = ["cmd", "/c", claude_cmd, "-p", "--output-format", "stream-json", "--verbose", "--model", "gpt-6-astra", "--max-turns", "200",
            "--permission-mode", "dontAsk", "--allowedTools", "Bash,PowerShell,Read,Write,Edit,Glob,Grep,NotebookEdit",
            "--disallowedTools", "WebFetch,WebSearch,Task,Agent", "--settings", str(CLAUDE_HOME / "settings.json"), "--setting-sources", "user"]
    t_start = time.monotonic()
    deadline = t_start + limit_s
    session_id = None
    last_text = ""
    replies_sent = {"no_info": 0, "scope": 0}
    attempts = []

    def one_attempt(argv, stdin_text, tag):
        """The prompt travels on stdin: a multi-line prompt on the command line is cut at the first newline by cmd."""
        nonlocal session_id, last_text
        rec = {"tag": tag, "argv": argv, "cwd": str(arena), "started_at": now(), "stdin_chars": len(stdin_text)}
        with stream.open("a", encoding="utf-8") as ev, err_path.open("a", encoding="utf-8") as err:
            proc = subprocess.Popen(argv, cwd=str(arena), env=env, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=err,
                                    text=True, encoding="utf-8", errors="replace", bufsize=1)
            ev.write(json.dumps({"type": "supervisor_start", "tag": tag, "pid": proc.pid, "ts": time.time()}) + "\n"); ev.flush()
            try:
                proc.stdin.write(stdin_text); proc.stdin.close()
            except OSError:
                pass
            reason = "completed"

            def pump():
                nonlocal session_id, last_text
                for line in proc.stdout:
                    ev.write(line if line.endswith("\n") else line + "\n"); ev.flush()
                    try:
                        e = json.loads(line)
                    except ValueError:
                        continue
                    if e.get("type") == "system" and e.get("session_id"):
                        session_id = e["session_id"]
                    if e.get("type") == "assistant":
                        for b in (e.get("message") or {}).get("content") or []:
                            if isinstance(b, dict) and b.get("type") == "text" and b.get("text"):
                                last_text = b["text"]
                    if e.get("type") == "result" and e.get("result"):
                        last_text = str(e["result"])
            th = threading.Thread(target=pump, daemon=True); th.start()
            while proc.poll() is None:
                if time.monotonic() > deadline:
                    reason = "timeout"
                    ev.write(json.dumps({"type": "supervisor_timeout", "ts": time.time()}) + "\n"); ev.flush()
                    kill_tree(proc.pid, lambda s: ev.write(s + "\n"))
                    break
                mk = stop_marker(run["run_id"])
                if mk:
                    reason = "stopped_by_controller"
                    early_stop.update(mk)
                    ev.write(json.dumps({"type": "supervisor_stop", "ts": time.time(), "marker": mk}, ensure_ascii=False) + "\n"); ev.flush()
                    kill_tree(proc.pid, lambda s: ev.write(s + "\n"))
                    break
                time.sleep(2)
            th.join(timeout=30)
            kill_tree(proc.pid, lambda s: ev.write(s + "\n"))
            rec.update({"ended_at": now(), "exit_code": proc.poll(), "exit_reason": reason})
        attempts.append(rec)
        return reason

    reason = one_attempt(base, ctx["prompt"], "initial")
    for _ in range(2):
        if reason != "completed" or time.monotonic() > deadline - 600 or not session_id or not looks_like_question(last_text):
            break
        if run["case_id"] == "nu1000" and mentions_scope(last_text) and replies_sent["scope"] == 0:
            reply, key = REPLY_SCOPE, "scope"
        elif replies_sent["no_info"] == 0:
            reply, key = REPLY_NO_INFO, "no_info"
        else:
            break
        replies_sent[key] += 1
        iv.add("fixed_reply", template=key, text=reply, trigger=last_text[-300:])
        reason = one_attempt(base + ["--resume", session_id], reply, f"resume_{key}")
    (evid / "commands" / "claude_attempts.json").write_text(json.dumps(attempts, ensure_ascii=False, indent=2), encoding="utf-8")
    (evid / "final_message.txt").write_text(last_text, encoding="utf-8")
    return {"exit_reason": reason, "session_id": session_id, "science_seconds": round(time.monotonic() - t_start, 1),
            "warmup_seconds": 0, "interventions_auto": iv.count_auto, "final_text": last_text[-2000:], "early_stop": early_stop or None}


def run_workbench(ctx: dict, run: dict, limit_s: int, knowledge_mode: str) -> dict:
    arena, evid = ctx["arena"], ctx["evid"]
    run_id = run["run_id"]
    iv = Interventions(evid / "interventions.jsonl")
    sink = Sink(evid / "raw_events" / "workbench_events.jsonl")
    wb = Workbench()
    proj = str(arena)
    # knowledge snapshot reset (H1 may have written skills in an earlier run)
    reset_knowledge(run_id)
    t_open = time.monotonic()
    opened = wb.open_project(proj)
    wb.set_settings(proj, permission_mode="auto", settings={**CONDITION_SETTINGS, "knowledge_mode": knowledge_mode})
    settings = wb.get_settings(proj)
    (evid / "effective_instructions" / "project_settings.json").write_text(json.dumps({"open": opened, "settings": settings}, ensure_ascii=False, indent=2), encoding="utf-8")
    if (arena / "AGENTS.md").exists():
        shutil.copy2(arena / "AGENTS.md", evid / "effective_instructions" / "AGENTS.md")
    # warm-up (non-scientific), then MCP readiness
    s = wb.send(proj, WARMUP, title=f"{run_id}")
    thread = s["thread_id"]
    st = wb.follow(thread, proj, sink, deadline=time.monotonic() + 900)
    mcp = wb.mcp_status(proj)
    if not (mcp.get("present") and (mcp.get("n_tools") or 0) > 0):
        log_event(run_id, "mcp_not_ready_retry", mcp=mcp)
        iv.add("engine_restart_for_mcp", detail=mcp)
        wb.restart_engine(proj)
        time.sleep(20)
        s = wb.send(proj, WARMUP, thread_id=thread)
        st = wb.follow(thread, proj, sink, cursor=st["cursor"], deadline=time.monotonic() + 900)
        mcp = wb.mcp_status(proj)
    warm = {"thread": thread, "ended": st["ended"], "mcp": mcp, "seconds": round(time.monotonic() - t_open, 1)}
    (evid / "commands" / "warmup.json").write_text(json.dumps(warm, ensure_ascii=False, indent=2), encoding="utf-8")
    if not (mcp.get("present") and (mcp.get("n_tools") or 0) > 0):
        sink.close()
        return {"exit_reason": "preflight_blocked", "detail": "MCP toolset not ready after one engine restart", "thread_id": thread,
                "warmup_seconds": warm["seconds"], "science_seconds": 0, "interventions_auto": iv.count_auto, "final_text": ""}

    # the science task: the clock starts here
    atts = [{"rel": "inputs/synthesis.png", "kind": "image", "name": "synthesis.png"}] if run["case_id"] == "nu1000" else None
    t_start = time.monotonic()
    deadline = t_start + limit_s
    state = {"last_text": "", "approvals_seen": set(), "bg_notice_sent": False, "bg_completed": False}
    replies_sent = {"no_info": 0, "scope": 0}

    def on_event(e: dict) -> None:
        k = e.get("kind")
        if k == "agent_message":
            state["last_text"] = e.get("text") or e.get("content") or state["last_text"]
        elif k == "approval_request":
            for a in wb.approvals(proj):
                aid = a.get("approval_id") or a.get("id")
                if aid and aid not in state["approvals_seen"]:
                    state["approvals_seen"].add(aid)
                    res = wb.decide(aid, "reject")
                    iv.add("approval_rejected", approval=a, result=res)
        elif k == "background_turn_completed":
            state["bg_completed"] = True

    wb.send(proj, ctx["prompt"], thread_id=thread, attachments=atts)
    cursor = st["cursor"]
    exit_reason = "completed"
    while True:
        st = wb.follow(thread, proj, sink, cursor=cursor, deadline=deadline, on_event=on_event)
        cursor = st["cursor"]
        if st["ended"] == "deadline":
            exit_reason = "timeout"
            res = wb.interrupt(thread)
            iv.add("interrupt_at_deadline", result=res)
            st = wb.follow(thread, proj, sink, cursor=cursor, deadline=time.monotonic() + 240, on_event=on_event)
            cursor = st["cursor"]
            break
        if st["ended"] == "turn_failed":
            exit_reason = "turn_failed"
            break
        if st["ended"] == "channel_closed":
            exit_reason = "channel_closed"
            break
        # turn completed: is the agent asking for facts it cannot have?
        text = state["last_text"]
        remaining = deadline - time.monotonic()
        if remaining > 600 and looks_like_question(text):
            if run["case_id"] == "nu1000" and mentions_scope(text) and replies_sent["scope"] == 0:
                reply, key = REPLY_SCOPE, "scope"
            elif replies_sent["no_info"] == 0:
                reply, key = REPLY_NO_INFO, "no_info"
            else:
                break
            replies_sent[key] += 1
            iv.add("fixed_reply", template=key, text=reply, trigger=text[-300:])
            wb.send(proj, reply, thread_id=thread)
            continue
        # a background job (e.g. SHELXT) may still be running after the agent's turn ended. The thread's busy
        # flag can lag the turn_completed event by a moment, so re-check a few times before committing to a wait
        # (rerun stayed in this wait for 30 min on a transient busy=True after the agent had already finished).
        def still_busy() -> bool:
            return any(t.get("busy") for t in wb.threads(proj).get("threads", [])) or has_running_jobs(arena)
        busy = False
        for _ in range(4):
            busy = still_busy()
            if not busy:
                break
            time.sleep(5)
        if remaining > 300 and busy and not state["bg_notice_sent"]:
            wait_until = min(deadline, time.monotonic() + 1800)
            st2 = {"cursor": cursor, "ended": "deadline"}
            while time.monotonic() < wait_until:
                st2 = wb.follow(thread, proj, sink, cursor=st2["cursor"], deadline=min(wait_until, time.monotonic() + 60), on_event=on_event,
                                stop_when=lambda e: e.get("kind") in ("background_turn_completed", "turn_started"))
                if st2["ended"] != "deadline" or not still_busy():
                    break
            cursor = st2["cursor"]
            if st2["ended"] in ("background_turn_completed", "turn_started") and remaining > 300:
                # the product itself surfaces the job result to the agent; only if it stays idle do we send the fixed notice
                if st2["ended"] == "turn_started":
                    continue
                state["bg_notice_sent"] = True
                iv.add("background_job_notice", text=NOTICE_BG)
                wb.send(proj, NOTICE_BG, thread_id=thread)
                continue
        break
    science_s = round(time.monotonic() - t_start, 1)
    # native transcript, all pages
    pages, before = [], None
    for _ in range(50):
        tr = wb.transcript(thread, proj, limit=2000, before=before)
        evs = tr.get("events", [])
        pages.append(tr)
        if not evs or len(evs) < 2000:
            break
        before = min(int(e.get("eid", 0)) for e in evs if isinstance(e, dict) and e.get("eid")) or None
        if not before:
            break
    (evid / "native_state" / "transcript_pages.json").write_text(json.dumps(pages, ensure_ascii=False), encoding="utf-8")
    try:
        (evid / "native_state" / "artifacts_index.json").write_text(json.dumps(wb.artifacts(thread, proj), ensure_ascii=False, indent=1), encoding="utf-8")
    except Exception as e:  # noqa: BLE001
        log_event(run_id, "artifacts_index_failed", error=str(e)[:200])
    sink.close()
    (evid / "final_message.txt").write_text(state["last_text"], encoding="utf-8")
    return {"exit_reason": exit_reason, "thread_id": thread, "warmup_seconds": warm["seconds"], "science_seconds": science_s,
            "interventions_auto": iv.count_auto, "final_text": state["last_text"][-2000:], "approvals_rejected": len(state["approvals_seen"])}


def has_running_jobs(arena: Path) -> bool:
    jobs = arena / ".crystalpilot" / "refine" / "jobs"
    if not jobs.is_dir():
        return False
    for j in jobs.glob("*/job.json"):
        try:
            d = json.loads(j.read_text(encoding="utf-8"))
            if str(d.get("status", "")).lower() in ("running", "started", "queued"):
                return True
        except Exception:  # noqa: BLE001
            continue
    return False


def reset_knowledge(run_id: str) -> None:
    """Restore the frozen H1 knowledge snapshot; archive anything a previous run wrote."""
    snap = json.loads((FREEZE / "knowledge_snapshot_manifest.json").read_text(encoding="utf-8"))
    current = manifest_dir(KNOWLEDGE)
    changed = [r for r in current if r not in snap or snap[r]["sha256"] != current[r]["sha256"]]
    removed = [r for r in snap if r not in current]
    if changed or removed:
        arch = CTRL / "provenance" / "knowledge_writes" / f"before_{run_id}"
        arch.mkdir(parents=True, exist_ok=True)
        for r in changed:
            dst = arch / r
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(KNOWLEDGE / r, dst)
            if r not in snap:
                (KNOWLEDGE / r).unlink()
        for r in list(changed) + removed:
            src = FREEZE / "knowledge_snapshot" / r
            if src.exists():
                (KNOWLEDGE / r).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, KNOWLEDGE / r)
        log_event(run_id, "knowledge_reset", changed=changed, removed=removed)


# ----------------------------------------------------------------------------- collection
def collect_and_seal(ctx: dict, run: dict, result: dict, limit_s: int) -> dict:
    arena, evid, run_id = ctx["arena"], ctx["evid"], run["run_id"]
    cond = run["condition"]
    # deliverables: the agent's declared delivery + the product results folder
    for name in ("deliverables", "CrystalPilot Results"):
        src = arena / name
        if src.exists():
            shutil.copytree(src, evid / "artifacts" / name, dirs_exist_ok=True)
    for name in ("AGENTS.md", "TASK.md", "environment.json", "SUMMARY.md", "README.md", "NOTES.md", "notes.md"):
        if (arena / name).exists():
            shutil.copy2(arena / name, evid / "artifacts" / name)
    if cond in ("H0", "H1"):
        cp = arena / ".crystalpilot"
        if (cp / "mcp_server.jsonl").exists():
            shutil.copy2(cp / "mcp_server.jsonl", evid / "native_state" / "mcp_server.jsonl")
        nodes = cp / "refine" / "nodes"
        if nodes.is_dir():
            for nj in nodes.glob("*/node.json"):
                dst = evid / "native_state" / "nodes" / nj.parent.name
                dst.mkdir(parents=True, exist_ok=True)
                shutil.copy2(nj, dst / "node.json")
        # anything the run wrote into the shared knowledge copy
        snap = json.loads((FREEZE / "knowledge_snapshot_manifest.json").read_text(encoding="utf-8"))
        current = manifest_dir(KNOWLEDGE)
        writes = [r for r in current if r not in snap or snap[r]["sha256"] != current[r]["sha256"]]
        (evid / "native_state" / "knowledge_writes.json").write_text(json.dumps(writes, ensure_ascii=False, indent=1), encoding="utf-8")
    # inputs unchanged?
    staged = ctx["manifest"]["inputs_manifest"]
    tampered = []
    for rel, meta in staged.items():
        p = arena / "inputs" / rel
        if not p.exists() or sha256_of(p) != meta["sha256"]:
            tampered.append(rel)
    seal = {"sealed_at": now(), "arena_manifest": manifest_dir(arena), "inputs_tampered": tampered,
            "evidence_manifest": manifest_dir(evid)}
    (evid / "seal_manifest.json").write_text(json.dumps(seal, ensure_ascii=False, indent=1), encoding="utf-8")
    usage = summarize_usage(evid, cond)
    (evid / "usage.json").write_text(json.dumps(usage, ensure_ascii=False, indent=2), encoding="utf-8")
    m = ctx["manifest"]
    m.update({"exit_reason": result.get("exit_reason"), "started_science_at": result.get("started_science_at"), "ended_at": now(),
              "warmup_seconds": result.get("warmup_seconds"), "science_seconds": result.get("science_seconds"),
              "interventions_auto": result.get("interventions_auto", 0), "human_scientific_interventions": 0,
              "inputs_tampered": tampered, "thread_or_session": result.get("thread_id") or result.get("session_id"),
              "early_stop": result.get("early_stop"),
              "protocol_validity": "valid" if not tampered else "protocol_invalid: inputs changed"})
    (evid / "manifest.json").write_text(json.dumps(m, ensure_ascii=False, indent=2), encoding="utf-8")
    return {"tampered": tampered, "usage": usage}


def summarize_usage(evid: Path, cond: str) -> dict:
    u = {"model_requests": None, "tool_calls": None, "input_tokens": None, "output_tokens": None, "cached_input_tokens": None,
         "reasoning_output_tokens": None, "cost_usd": None, "source": None}
    if cond in ("C0", "C0B", "C0N"):
        p = evid / "raw_events" / "codex_events.jsonl"
        tools = 0; inp = out = cached = reas = 0; turns = 0; msgs = 0
        for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
            try:
                e = json.loads(line)
            except ValueError:
                continue
            it = e.get("item") or {}
            if e.get("type") == "item.completed" and it.get("type") in ("command_execution", "mcp_tool_call", "file_change", "web_search"):
                tools += 1
            if e.get("type") == "item.completed" and it.get("type") in ("agent_message", "reasoning"):
                msgs += 1
            if e.get("type") == "turn.completed":
                turns += 1
                us = e.get("usage") or {}
                inp += us.get("input_tokens", 0); out += us.get("output_tokens", 0)
                cached += us.get("cached_input_tokens", 0); reas += us.get("reasoning_output_tokens", 0)
        u.update({"tool_calls": tools, "input_tokens": inp, "output_tokens": out, "cached_input_tokens": cached, "reasoning_output_tokens": reas,
                  "turns": turns, "source": "codex exec JSONL: turn.completed usage totals; model_requests not exposed (null)"})
    elif cond in ("P1", "P1B", "P1N"):
        p = evid / "raw_events" / "claude_stream.jsonl"
        tools = 0; req = 0; res = None
        for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
            try:
                e = json.loads(line)
            except ValueError:
                continue
            if e.get("type") == "assistant":
                req += 1
                for b in (e.get("message") or {}).get("content") or []:
                    if isinstance(b, dict) and b.get("type") == "tool_use":
                        tools += 1
            if e.get("type") == "result":
                res = e
        us = (res or {}).get("usage") or {}
        u.update({"model_requests": req, "tool_calls": tools, "input_tokens": us.get("input_tokens"), "output_tokens": us.get("output_tokens"),
                  "cached_input_tokens": us.get("cache_read_input_tokens"), "cost_usd": (res or {}).get("total_cost_usd"),
                  "num_turns": (res or {}).get("num_turns"), "source": "claude stream-json: assistant messages counted as requests; result.usage totals; cost is Claude Code's own estimate for an unknown model and is not reliable"})
    else:
        p = evid / "raw_events" / "workbench_events.jsonl"
        tools = 0; inp = out = cached = reas = 0; n_usage = 0; last = None
        for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
            try:
                e = json.loads(line).get("event") or {}
            except ValueError:
                continue
            k = e.get("kind") or ""
            # every tool-like item the workbench surfaces starts with a *_started event (shell command,
            # file change, MCP tool); turn_started is the turn itself
            if k.endswith("_started") and k != "turn_started":
                tools += 1
            if k == "token_usage":
                n_usage += 1
                last = e
        if last:
            tu = last.get("total") or last.get("usage") or last
            inp = tu.get("input_tokens") or tu.get("inputTokens"); out = tu.get("output_tokens") or tu.get("outputTokens")
            cached = tu.get("cached_input_tokens") or tu.get("cachedInputTokens"); reas = tu.get("reasoning_output_tokens") or tu.get("reasoningOutputTokens")
        u.update({"tool_calls": tools, "input_tokens": inp, "output_tokens": out, "cached_input_tokens": cached, "reasoning_output_tokens": reas,
                  "token_usage_events": n_usage, "last_token_usage_event": last, "source": "workbench token_usage events (last cumulative record) and tool_started count; model_requests not exposed (null)"})
    return u


# ----------------------------------------------------------------------------- main
def main() -> int:
    cap_note = cpucap.limit_self()
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--limit-minutes", type=int, default=None, help="mechanics smoke test only")
    ap.add_argument("--plan", default=str(FREEZE / "plan.json"))
    ap.add_argument("--allow-concurrent", action="store_true",
                    help="user-requested extra trial that runs beside the driver's current run; recorded in the manifest")
    a = ap.parse_args()
    plan = json.loads(Path(a.plan).read_text(encoding="utf-8"))
    run = next((r for r in plan["runs"] if r["run_id"] == a.run_id), None)
    if run is None:
        raise SystemExit(f"unknown run id {a.run_id}")
    limit_s = (a.limit_minutes or plan["cases"][run["case_id"]]["limit_minutes"]) * 60
    concurrent_with = None
    lock = LOCK
    if LOCK.exists():
        if not a.allow_concurrent:
            raise SystemExit(f"another run is active: {LOCK.read_text()}")
        concurrent_with = LOCK.read_text(encoding="utf-8").strip()
        lock = LOCK.with_name(f"active_run_{a.run_id}.lock")
    free_gb = shutil.disk_usage("H:\\").free / 1e9
    if free_gb < 5:
        raise SystemExit(f"disk too low: {free_gb:.1f} GB")
    lock.write_text(a.run_id, encoding="utf-8")
    try:
        log_event(a.run_id, "preparing", condition=run["condition"], case=run["case_id"], limit_minutes=limit_s // 60, disk_free_gb=round(free_gb, 1), cpu_cap=cap_note,
                  concurrent_with=concurrent_with, plan=Path(a.plan).name)
        update_state(a.run_id, status="preparing", condition=run["condition"], case_id=run["case_id"], limit_seconds=limit_s, started_at=now(),
                     concurrent_with=concurrent_with, plan=Path(a.plan).name, extra=run.get("extra"))
        ctx = prepare(run, plan, limit_s)
        ctx["manifest"]["concurrent_with"] = concurrent_with
        ctx["manifest"]["plan_file"] = Path(a.plan).name
        ctx["manifest"]["extra"] = run.get("extra")
        (ctx["evid"] / "manifest.json").write_text(json.dumps(ctx["manifest"], ensure_ascii=False, indent=2), encoding="utf-8")
        update_state(a.run_id, status="running")
        log_event(a.run_id, "started", arena=str(ctx["arena"]))
        t0 = now()
        if run["condition"] in ("C0", "C0B", "C0N"):
            result = run_c0(ctx, run, limit_s, bare=run["condition"] == "C0B", pure=run["condition"] == "C0N")
            result["orphans_killed"] = cpucap.kill_run_orphans(ctx["arena"])
        elif run["condition"] in ("P1", "P1B", "P1N"):
            result = run_p1(ctx, run, limit_s, bare=run["condition"] == "P1B", pure=run["condition"] == "P1N")
            result["orphans_killed"] = cpucap.kill_run_orphans(ctx["arena"])
        else:
            result = run_workbench(ctx, run, limit_s, "tools_only" if run["condition"] == "H0" else "full")
        result["started_science_at"] = t0
        log_event(a.run_id, "backend_finished", exit_reason=result.get("exit_reason"), science_seconds=result.get("science_seconds"))
        update_state(a.run_id, status="collecting", exit_reason=result.get("exit_reason"))
        sealed = collect_and_seal(ctx, run, result, limit_s)
        status = "sealed" if not sealed["tampered"] else "sealed_protocol_invalid"
        update_state(a.run_id, status=status, ended_at=now(), science_seconds=result.get("science_seconds"),
                     warmup_seconds=result.get("warmup_seconds"), interventions_auto=result.get("interventions_auto"),
                     usage=sealed["usage"], evidence=str(ctx["evid"]))
        log_event(a.run_id, "sealed", status=status, usage={k: sealed["usage"].get(k) for k in ("tool_calls", "input_tokens", "output_tokens")})
        print("FINAL TEXT (tail):", (result.get("final_text") or "")[-800:])
        return 0
    except Exception as e:  # noqa: BLE001
        log_event(a.run_id, "controller_error", error=f"{type(e).__name__}: {str(e)[:500]}")
        update_state(a.run_id, status="infra_failed", error=f"{type(e).__name__}: {str(e)[:500]}", ended_at=now())
        raise
    finally:
        lock.unlink(missing_ok=True)


if __name__ == "__main__":
    sys.exit(main())
