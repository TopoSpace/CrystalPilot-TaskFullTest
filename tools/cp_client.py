"""Thin REST/SSE client for the experiment's CrystalPilot workbench instance (port 8021).

Every server event is appended to a JSONL sink as received (seq, wall time, payload), so the trial's
native stream survives a controller crash. Reconnects resume from the last seq; the channel's
`generation` is recorded so a server restart can be recognised and the transcript route used to fill gaps.
"""
import json
import time
from pathlib import Path
from typing import Callable, Optional

import httpx

BASE = "http://127.0.0.1:8021/api"
TURN_END = ("turn_completed", "turn_failed", "turn_aborted", "turn_interrupted")


class Sink:
    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path = path
        self.fh = path.open("a", encoding="utf-8")

    def write(self, seq: int, payload: dict, note: str | None = None) -> None:
        rec = {"seq": seq, "wall": time.strftime("%Y-%m-%dT%H:%M:%S"), "ts": round(time.time(), 3), "event": payload}
        if note:
            rec["note"] = note
        self.fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
        self.fh.flush()

    def close(self) -> None:
        self.fh.close()


class Workbench:
    def __init__(self, base: str = BASE, timeout: float = 120.0):
        self.base = base
        self.http = httpx.Client(timeout=httpx.Timeout(timeout, connect=15))

    # -- projects -----------------------------------------------------------
    def health(self) -> dict:
        return self.http.get(self.base + "/health").json()

    def open_project(self, path: str) -> dict:
        r = self.http.post(self.base + "/projects/open", json={"path": path, "auto_approve": False})
        r.raise_for_status()
        return r.json()

    def get_settings(self, path: str) -> dict:
        r = self.http.get(self.base + "/projects/settings", params={"path": path})
        r.raise_for_status()
        return r.json()

    def set_settings(self, path: str, permission_mode: str | None = None, settings: dict | None = None) -> dict:
        body: dict = {"path": path}
        if permission_mode:
            body["permission_mode"] = permission_mode
        if settings is not None:
            body["settings"] = settings
        r = self.http.post(self.base + "/projects/settings", json=body)
        r.raise_for_status()
        return r.json()

    def mcp_status(self, path: str) -> dict:
        r = self.http.post(self.base + "/projects/mcp_status", json={"path": path}, timeout=300)
        r.raise_for_status()
        return r.json()

    def restart_engine(self, path: str) -> dict:
        r = self.http.post(self.base + "/projects/restart_engine", json={"path": path}, timeout=300)
        r.raise_for_status()
        return r.json()

    def usage(self, path: str) -> dict:
        r = self.http.get(self.base + "/projects/usage", params={"path": path})
        return r.json()

    # -- threads ------------------------------------------------------------
    def send(self, project: str, message: str, thread_id: str | None = None, title: str | None = None,
             attachments: list | None = None) -> dict:
        body: dict = {"project": project, "message": message}
        if thread_id:
            body["thread_id"] = thread_id
        if title:
            body["title"] = title
        if attachments:
            body["attachments"] = attachments
        r = self.http.post(self.base + "/threads/send", json=body)
        r.raise_for_status()
        return r.json()

    def interrupt(self, thread_id: str) -> dict:
        r = self.http.post(self.base + "/threads/interrupt", json={"thread_id": thread_id}, timeout=60)
        return {"status": r.status_code, "body": r.text[:500]}

    def threads(self, project: str) -> dict:
        return self.http.get(self.base + "/threads/list", params={"project": project}).json()

    def transcript(self, thread_id: str, project: str, limit: int = 2000, before: int | None = None) -> dict:
        params: dict = {"thread_id": thread_id, "project": project, "limit": limit}
        if before is not None:
            params["before"] = before
        return self.http.get(self.base + "/threads/transcript", params=params, timeout=300).json()

    def approvals(self, path: str) -> list:
        r = self.http.get(self.base + "/approvals", params={"path": path})
        return r.json() if r.status_code == 200 else []

    def decide(self, approval_id: str, decision: str) -> dict:
        r = self.http.post(self.base + "/approvals/decide", json={"approval_id": approval_id, "decision": decision})
        return {"status": r.status_code, "body": r.text[:300]}

    def artifacts(self, thread_id: str, project: str) -> dict:
        r = self.http.get(self.base + "/threads/artifacts", params={"thread_id": thread_id, "project": project})
        return r.json() if r.status_code == 200 else {"status": r.status_code}

    # -- event stream ---------------------------------------------------------
    def follow(self, thread_id: str, project: str, sink: Sink, cursor: int = 0, *,
               deadline: float, on_event: Optional[Callable[[dict], None]] = None,
               stop_when: Optional[Callable[[dict], bool]] = None) -> dict:
        """Stream events into `sink` until a turn-end event (or stop_when) fires or `deadline`
        (time.monotonic) passes. Returns cursor, last kind, generation and how it ended."""
        state = {"cursor": cursor, "last_kind": None, "generation": None, "ended": None, "reconnects": 0}
        while time.monotonic() < deadline:
            try:
                with self.http.stream("GET", self.base + "/threads/events",
                                      params={"thread_id": thread_id, "after": state["cursor"], "project": project},
                                      timeout=httpx.Timeout(75, connect=15)) as resp:
                    resp.raise_for_status()
                    seq, data = state["cursor"], []
                    for line in resp.iter_lines():
                        if line.startswith("id:"):
                            seq = int(line[3:].strip())
                        elif line.startswith("data:"):
                            data.append(line[5:].lstrip())
                        elif not line and data:
                            payload = json.loads("\n".join(data))
                            data = []
                            kind = payload.get("kind")
                            if kind == "channel_hello":
                                gen = payload.get("generation")
                                if state["generation"] is not None and gen != state["generation"]:
                                    sink.write(seq, payload, note="generation changed: server restarted; fill from transcript")
                                state["generation"] = gen
                                continue
                            if kind == "ping":
                                continue
                            if seq <= state["cursor"] and kind != "channel_closed":
                                continue
                            state["cursor"] = max(state["cursor"], seq)
                            state["last_kind"] = kind
                            sink.write(seq, payload)
                            if on_event:
                                on_event(payload)
                            if kind in TURN_END or (stop_when and stop_when(payload)):
                                state["ended"] = kind
                                return state
                            if kind == "channel_closed":
                                state["ended"] = "channel_closed"
                                return state
                        if time.monotonic() >= deadline:
                            break
            except (httpx.HTTPError, ValueError) as err:
                state["reconnects"] += 1
                sink.write(state["cursor"], {"kind": "client_reconnect", "error": type(err).__name__, "detail": str(err)[:200]})
                time.sleep(3)
        state["ended"] = "deadline"
        return state
