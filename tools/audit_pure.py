"""Rule audit of the pure-baseline runs (C0N / P1N): every shell command and every file-tool call is checked for
crystallographic programs or libraries, other interpreters or environments, downloads or installs, and for paths
outside the run's own arena (the plain Python and the temp directory excepted). Read-only; writes
control/analysis/pure_audit.json and prints a summary."""
import json
import re
import sys
from pathlib import Path

ROOT = Path(r"H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921")
CTRL = ROOT / "control"
FORBIDDEN = re.compile(r"shelx|olex2|platon|dials|\bxds\b|crysalis|cctbx|gemmi|pymatgen|\base\b|spglib|crystalpilot(?!-campaigns)|"
                       r"miniforge|pyenv|conda|pip install|pip3 install|winget|npm install|Invoke-WebRequest|curl |wget |git clone|Start-BitsTransfer|"
                       r"D:\\OLEX2|H:\\CrystalPilot\\", re.I)
#: the plain Python, temp dirs, Windows system binaries (powershell.exe itself) and the isolated Claude Code home
#: (its own background-task output files live there) are part of the run's environment, not outside access
ALLOWED_PREFIXES = ("toolbox\\pyplain", "c:\\tmp\\", "c:\\windows\\", "engine\\claude-home\\")


def collect(evid: Path) -> tuple:
    cmds, files = [], []
    cs, cx = evid / "raw_events" / "claude_stream.jsonl", evid / "raw_events" / "codex_events.jsonl"
    if cs.exists():
        for l in cs.read_text(encoding="utf-8", errors="replace").splitlines():
            try:
                e = json.loads(l)
            except ValueError:
                continue
            if e.get("type") != "assistant":
                continue
            for b in (e.get("message") or {}).get("content") or []:
                if not (isinstance(b, dict) and b.get("type") == "tool_use"):
                    continue
                inp = b.get("input") or {}
                if b.get("name") in ("Read", "Glob", "Grep", "Edit", "Write", "NotebookEdit"):
                    files.append((b["name"], str(inp.get("file_path") or inp.get("path") or inp.get("pattern") or "")))
                else:
                    files_in = str(inp.get("command") or json.dumps(inp, ensure_ascii=False))
                    cmds.append((b.get("name"), files_in))
    if cx.exists():
        for l in cx.read_text(encoding="utf-8", errors="replace").splitlines():
            try:
                e = json.loads(l)
            except ValueError:
                continue
            it = e.get("item") or {}
            if e.get("type") == "item.completed" and it.get("type") == "command_execution":
                cmds.append(("shell", it.get("command") or ""))
    return cmds, files


def outside_arena(path: str, arena: str) -> bool:
    p = path.lower().replace("/", "\\")
    while "\\\\" in p:                    # doubled backslashes from PowerShell quoting inside command strings
        p = p.replace("\\\\", "\\")
    if not re.match(r"[a-z]:\\", p):
        return False                      # relative path: inside the working directory
    if p.startswith(arena):
        return False
    return not any(a in p for a in ALLOWED_PREFIXES)


def main() -> int:
    rids = sys.argv[1:] or ["r09_P1N_nu1000", "r10_C0N_nu1000", "r11_C0N_alanine", "r12_P1N_alanine"]
    out = {}
    for rid in rids:
        evid = CTRL / "runs" / rid
        arena = str(ROOT / "arena" / rid).lower()
        cmds, files = collect(evid)
        hits = [(n, c[:200]) for n, c in cmds if FORBIDDEN.search(c)]
        abs_paths = []
        for n, c in cmds:
            for m in re.findall(r"[A-Za-z]:\\[^\s\"'`;)]+", c):
                if outside_arena(m, arena):
                    abs_paths.append((n, m[:160]))
        file_out = [(n, p[:160]) for n, p in files if outside_arena(p, arena)]
        interp = sorted({m.lower().replace("\\\\", "\\") for n, c in cmds for m in re.findall(r"[A-Za-z]:\\[^\s\"'`;)]*python[^\s\"'`;)]*\.exe", c, re.I)})
        rec = {"shell_commands": len(cmds), "file_tool_calls": len(files), "forbidden_pattern_hits": hits, "commands_with_paths_outside_arena": abs_paths,
               "file_tool_calls_outside_arena": file_out, "python_executables_referenced": interp,
               "verdict": "clean" if not hits and not abs_paths and not file_out and all("toolbox\\pyplain" in i for i in interp) else "review"}
        out[rid] = rec
        print(f"== {rid}: {rec['verdict']} | commands {len(cmds)} | file-tool calls {len(files)} | forbidden hits {len(hits)} | "
              f"outside-arena command paths {len(abs_paths)} | outside-arena file-tool paths {len(file_out)} | interpreters {interp}")
        for n, c in hits[:5]:
            print("   HIT:", n, c.replace("\n", " "))
        for n, c in abs_paths[:5]:
            print("   PATH:", n, c)
        for n, c in file_out[:5]:
            print("   FILE:", n, c)
    (CTRL / "analysis" / "pure_audit.json").write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
