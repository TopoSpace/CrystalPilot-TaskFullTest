"""Write control/environment/environment.json: versions, tool inventory, machine facts, hashes of the
engines every condition shares. Keys and gateway host are never written; the provider is described by id only."""
import hashlib, json, os, platform, shutil, subprocess, sys, time, tomllib
from pathlib import Path

ROOT = Path(r"H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921")
ENGINE = ROOT / "engine" / "CrystalPilot"
MAIN = Path(r"H:\CrystalPilot")
DIALS = Path(r"C:\users\<user>\miniforge3\envs\dials")
OUT = ROOT / "control" / "environment"

def sha(p): 
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""): h.update(c)
    return h.hexdigest()

def run(cmd, timeout=60, **kw):
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, encoding="utf-8", errors="replace", **kw)
        return (p.stdout + p.stderr).strip()
    except Exception as e:
        return f"ERR {type(e).__name__}: {e}"

git = {"commit": run(["git", "-C", str(MAIN), "rev-parse", "HEAD"]), "short": run(["git", "-C", str(MAIN), "rev-parse", "--short", "HEAD"]),
       "status_lines": len([l for l in run(["git", "-C", str(MAIN), "status", "--porcelain"]).splitlines() if l.strip()]),
       "head_subject": run(["git", "-C", str(MAIN), "log", "-1", "--format=%s"])}
codex_bin = ENGINE / "vendor/codex/node_modules/@openai/codex-win32-x64/vendor/x86_64-pc-windows-msvc/bin/codex.exe"
kernel = {"binary": str(codex_bin), "version": run([str(codex_bin), "--version"]), "sha256": sha(codex_bin),
          "exec_help_saved": "kernel_exec_help.txt"}
(OUT / "kernel_exec_help.txt").write_text(run([str(codex_bin), "exec", "--help"]), encoding="utf-8")
claude = {"path": shutil.which("claude"), "version": run(["claude", "--version"], shell=True), "help_saved": "claude_help.txt"}
(OUT / "claude_help.txt").write_text(run(["claude", "--help"], shell=True), encoding="utf-8")
venv_py = MAIN / ".venv/Scripts/python.exe"
pkgs = run([str(venv_py), "-c", "import importlib,json\nout={}\nfor m in ['cctbx','iotbx','smtbx','gemmi','numpy','scipy','fastapi','uvicorn','httpx','psutil','pydantic']:\n  try: out[m]=getattr(importlib.import_module(m),'__version__','present')\n  except Exception as e: out[m]='missing'\nprint(json.dumps(out))"])
dials_pkgs = run([str(DIALS / "python.exe"), "-c", "import importlib,json,sys\nout={'python':sys.version.split()[0]}\nfor m in ['dials','dxtbx','cctbx','iotbx','smtbx','gemmi','numpy','scipy','matplotlib','PIL']:\n  try: out[m]=getattr(importlib.import_module(m),'__version__','present')\n  except Exception: out[m]='missing'\nprint(json.dumps(out))"], cwd=str(ROOT))
shelx = ENGINE / "vendor/shelx"
tools = {
    "shelxl": {"path": str(shelx / "shelxl.exe"), "sha256": sha(shelx / "shelxl.exe"), "banner": "SHELXL 2019/3 (multi-CPU), version banner captured on 2026-09-21"},
    "shelxt": {"path": str(shelx / "shelxt.exe"), "sha256": sha(shelx / "shelxt.exe"), "banner": "SHELXT 2018/2"},
    "platon": {"path": str(shelx / "platon.exe"), "sha256": sha(shelx / "platon.exe"), "note": "PLATON build dated 2024-01-11 by file time; banner not captured (interactive prompt)"},
    "olex2": {"path": str(ENGINE / "vendor/olex2/app/olex2.exe"), "exists": (ENGINE / "vendor/olex2/app/olex2.exe").exists()},
    "dials": {"prefix": str(DIALS), "version": run([str(DIALS / "Scripts/dials.version.exe")], env={**os.environ, "PATH": f"{DIALS/'Library/bin'};{DIALS/'Scripts'};{DIALS};{os.environ['PATH']}"}).splitlines()[:1],
              "dispatchers": str(DIALS / "Scripts"), "dll_dir": str(DIALS / "Library/bin"),
              "note": "the DIALS conda environment (Python 3.14, cctbx 2026.7) is also the plain scientific Python offered to the baseline conditions; it cannot import the CrystalPilot package"},
    "crystalpilot_venv": {"python": str(venv_py), "packages": json.loads(pkgs) if pkgs.startswith("{") else pkgs,
                          "note": "used only by the CrystalPilot server and MCP tools (H0/H1); the crystalpilot package is pip-installed editable in it, so it is not offered to C0/P1"},
    "dials_env_packages": json.loads(dials_pkgs) if dials_pkgs.startswith("{") else dials_pkgs,
}
cfg = tomllib.loads((ROOT / "engine/codex-home-wb/config.toml").read_text(encoding="utf-8"))
model = {"model_id": cfg["model"], "provider_id": cfg["model_provider"], "reasoning_effort": cfg["model_reasoning_effort"],
         "wire_api": cfg["model_providers"]["crystalpilot"]["wire_api"], "auth": "auth.command hook (bearer token supplied on demand by the product's credential script; key file outside the experiment tree)",
         "catalog_entry": next((m for m in json.load(open(ROOT / "engine/codex-home-wb/model_catalog.custom.json", encoding="utf-8"))["models"] if m.get("slug") == cfg["model"]), None),
         "gateway_probe_2026-09-21": {"responses_endpoint": "HTTP 200 for gpt-6-astra on both /openai and /openai/v1 base paths",
                                      "models_endpoint": "HTTP 200; gpt-6-astra listed with supported_endpoints /responses, reasoning_effort ladder low..max, vision, 1.178M context",
                                      "anthropic_messages_endpoint": "/openai/v1/messages answers HTTP 200 and reports model=gpt-6-astra when no effort control is sent; sending output_config.effort produced an error that names gpt-5-mini, so effort cannot be steered on that route",
                                      "chat_completions": "gpt-6-astra not offered on /chat/completions"}}
mem = run(["powershell", "-NoProfile", "-Command", "(Get-CimInstance Win32_OperatingSystem).TotalVisibleMemorySize"])
machine = {"os": platform.platform(), "logical_cpus": os.cpu_count(), "ram_gb": round(int(mem) / 1048576, 1) if mem.isdigit() else mem,
           "disk_free_gb": {d: round(shutil.disk_usage(d).free / 1e9, 1) for d in ("H:\\", "C:\\")}, "timezone": time.strftime("%Z %z"),
           "cpu_cap": "CRYSTALPILOT_CPU_CORES=4 (BelowNormal priority) for the workbench engine; one science task at a time"}
env = {"recorded_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "git": git, "engine_copy": str(ENGINE), "kernel": kernel, "claude_code": claude,
       "model": model, "tools": tools, "machine": machine,
       "experiment_server": {"port": 8021, "codex_home": str(ROOT / "engine/codex-home-wb"), "preferences": "engine copy's own workdir/preferences.json"},
       "production_server_note": "port 8010 was not listening at 14:35 on 2026-09-21; it was left alone"}
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "environment.json").write_text(json.dumps(env, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps({k: env[k] for k in ("git", "kernel", "claude_code", "machine")}, ensure_ascii=False, indent=1)[:1500])
print("dials:", tools["dials"]["version"], "| venv pkgs:", tools["crystalpilot_venv"]["packages"])
