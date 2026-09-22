"""CPU cap for the experiment's own process trees (controller, runners, verifier, native kernels).

Windows children inherit the parent's affinity mask, so capping a tree root (driver, run_trial, verify_run)
caps every codex / claude / SHELXL / SHELXT / PLATON descendant. Default: the first CRYSTALPILOT_EXPERIMENT_CORES
cores (8 of 32 = 25 %; user limit 2026-09-21 23:45: never above 30 %) at BelowNormal priority. The experiment
workbench server keeps its own 4-core cap (cores 0-3), a subset of these cores, so the two trees together never
exceed the 8 cores.

Usage:
  import cpucap; cpucap.limit_self()          at the top of a tree root
  python cpucap.py <pid> [<pid> ...]          clamp running processes
  python cpucap.py --sweep                    clamp every process under the controller's recorded roots and
                                              kill orphaned experiment-toolbox processes (PLATON never exits
                                              by itself; a "platon -h" left by an agent's shell spins one core
                                              forever). Only processes whose executable lives in the experiment
                                              toolbox (or the vendor directory it points to) AND whose parent is
                                              gone are killed; never by name.
"""
from __future__ import annotations

import os
import re
import sys
import time
from pathlib import Path

ROOT = Path(r"H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921")
CTRL = ROOT / "control"
TOOLBOX = ROOT / "toolbox"
BELOW_NORMAL = 0x00004000
PID_FILES = ("driver.pid", "chain_after.pid", "driver_bare.pid", "extra_rerun.pid")


def cores() -> int:
    try:
        n = int(os.environ.get("CRYSTALPILOT_EXPERIMENT_CORES", "8"))
    except ValueError:
        n = 8
    return max(1, min(n, os.cpu_count() or n))


def _k32():
    import ctypes
    import ctypes.wintypes as wt
    k = ctypes.windll.kernel32
    k.GetCurrentProcess.restype = ctypes.c_void_p
    k.OpenProcess.argtypes = [wt.DWORD, wt.BOOL, wt.DWORD]
    k.OpenProcess.restype = ctypes.c_void_p
    k.CloseHandle.argtypes = [ctypes.c_void_p]
    k.SetPriorityClass.argtypes = [ctypes.c_void_p, wt.DWORD]
    k.SetPriorityClass.restype = wt.BOOL
    k.GetPriorityClass.argtypes = [ctypes.c_void_p]
    k.GetPriorityClass.restype = wt.DWORD
    k.SetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
    k.SetProcessAffinityMask.restype = wt.BOOL
    k.GetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t), ctypes.POINTER(ctypes.c_size_t)]
    k.GetProcessAffinityMask.restype = wt.BOOL
    return k


def state(pid: int) -> dict | None:
    import ctypes
    k = _k32()
    h = k.OpenProcess(0x1000, False, int(pid))
    if not h:
        return None
    try:
        pm, sm = ctypes.c_size_t(), ctypes.c_size_t()
        ok = k.GetProcessAffinityMask(h, ctypes.byref(pm), ctypes.byref(sm))
        return {"mask": int(pm.value) if ok else None, "priority": int(k.GetPriorityClass(h))}
    finally:
        k.CloseHandle(h)


def limit(pid: int | None = None, n: int | None = None) -> str:
    """Cap one process (None = the current one). Returns what the kernel reports afterwards."""
    if os.name != "nt":
        return "not windows; no cap"
    import ctypes
    n = n or cores()
    k = _k32()
    if pid is None:
        h, close = k.GetCurrentProcess(), False
    else:
        h = k.OpenProcess(0x0200 | 0x1000, False, int(pid))
        if not h:
            return f"pid {pid}: cannot open"
        close = True
    try:
        mask = (1 << n) - 1
        ok_p = bool(k.SetPriorityClass(h, BELOW_NORMAL))
        ok_a = bool(k.SetProcessAffinityMask(h, mask))
        pm, sm = ctypes.c_size_t(), ctypes.c_size_t()
        k.GetProcessAffinityMask(h, ctypes.byref(pm), ctypes.byref(sm))
        prio = int(k.GetPriorityClass(h))
    finally:
        if close:
            k.CloseHandle(h)
    who = f"pid {pid}" if pid else "self"
    if ok_p and ok_a and pm.value == mask and prio == BELOW_NORMAL:
        return f"{who}: cpu-limited to {n} core(s) at BelowNormal"
    return f"{who}: cap NOT fully applied (priority ok={ok_p}, affinity ok={ok_a}, mask now {pm.value:#x}, priority {prio:#x})"


def limit_self() -> str:
    return limit(None)


def _toolbox_dirs() -> list:
    dirs = []
    for sub in ("shelx", "olex2"):
        p = TOOLBOX / sub
        if p.exists():
            dirs.append(str(p).lower().rstrip("\\"))
            try:
                dirs.append(str(p.resolve()).lower().rstrip("\\"))   # the vendor directory behind the junction
            except OSError:
                pass
    return dirs


def recorded_roots() -> list:
    pids = []
    for name in PID_FILES:
        p = CTRL / "controller" / name
        if p.exists():
            m = re.search(r"pid=(\d+)", p.read_text(encoding="utf-8", errors="replace"))
            if m:
                pids.append(int(m.group(1)))
    return pids


def sweep(kill_orphans: bool = True) -> list:
    """Clamp every live process under the recorded roots; kill orphaned toolbox processes."""
    import psutil
    n = cores()
    mask = (1 << n) - 1
    notes = []
    seen = set()
    for root in recorded_roots():
        try:
            rp = psutil.Process(root)
        except psutil.Error:
            continue
        for p in [rp] + rp.children(recursive=True):
            if p.pid in seen:
                continue
            seen.add(p.pid)
            st = state(p.pid)
            if st and (st["mask"] != mask or st["priority"] != BELOW_NORMAL):
                try:
                    name = p.name()
                except psutil.Error:
                    name = "?"
                notes.append(f"clamp {name} " + limit(p.pid, n))
    if kill_orphans:
        dirs = _toolbox_dirs()
        for p in psutil.process_iter(["pid", "ppid", "exe", "create_time", "name"]):
            try:
                exe = (p.info.get("exe") or "").lower()
                if not exe or not any(exe.startswith(d + "\\") for d in dirs):
                    continue
                if time.time() - (p.info.get("create_time") or time.time()) < 120:
                    continue
                parent_alive = psutil.pid_exists(p.info["ppid"]) and p.info["ppid"] not in (0, 4)
                if parent_alive:
                    # parent still there: only clamp it (its owner may still be reading its output)
                    st = state(p.pid)
                    if st and st["mask"] != mask:
                        notes.append(f"clamp orphan-candidate {p.info['name']} " + limit(p.pid, n))
                    continue
                cmd = " ".join(p.cmdline())[:120]
                p.kill()
                notes.append(f"killed orphan {p.info['name']} pid {p.pid} (parent {p.info['ppid']} gone, started {time.strftime('%H:%M:%S', time.localtime(p.info['create_time']))}, cmd {cmd})")
            except (psutil.Error, OSError):
                continue
    return notes


def kill_run_orphans(arena: Path) -> list:
    """After a native-kernel run: kill experiment-toolbox processes whose working directory is inside the
    run's arena (PLATON left behind by the agent's shells). Nothing outside the arena is touched."""
    import psutil
    dirs = _toolbox_dirs()
    arena_l = str(arena).lower().rstrip("\\")
    notes = []
    for p in psutil.process_iter(["pid", "exe", "name"]):
        try:
            exe = (p.info.get("exe") or "").lower()
            if not exe or not any(exe.startswith(d + "\\") for d in dirs):
                continue
            cwd = (p.cwd() or "").lower().rstrip("\\")
            if cwd == arena_l or cwd.startswith(arena_l + "\\"):
                p.kill()
                notes.append(f"killed {p.info['name']} pid {p.pid} (cwd in arena)")
        except (psutil.Error, OSError):
            continue
    return notes


if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0] == "--sweep":
        out = sweep()
        print("\n".join(out) if out else "sweep: nothing to change")
    else:
        for arg in args:
            print(limit(int(arg)))
