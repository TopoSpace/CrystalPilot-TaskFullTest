"""Non-scientific readiness probe of the experiment workbench for one knowledge mode.

Opens a scratch project, applies the frozen condition settings, runs the warm-up turn the trials will
use, lists the MCP tools actually exposed (frozen as the condition's tool schema), and checks with a
synthetic picture that an image attachment really reaches the model. No scientific input is touched.
"""
import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from cp_client import Sink, Workbench  # noqa: E402

ROOT = Path(r"H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921")
WARMUP = "不开始任何任务，也不要读取任何文件。只回复一行“就绪”。"
IMAGE_Q = "请只描述这张图片里画了什么（形状、颜色、文字），一句话，不要做其他事情。"

CONDITION_SETTINGS = {
    "subagents": "off",
    "enable_specialists": False,
    "allow_iucr_upload": False,
    "model_override": "gpt-6-astra",
    "effort_override": "xhigh",
    "model_provider_override": "crystalpilot",
}


def make_picture(path: Path) -> None:
    from PIL import Image, ImageDraw
    img = Image.new("RGB", (320, 200), (250, 250, 250))
    d = ImageDraw.Draw(img)
    d.ellipse((30, 30, 170, 170), fill=(200, 30, 30))
    d.rectangle((200, 60, 300, 160), fill=(30, 60, 200))
    d.text((40, 180), "PROBE 7", fill=(0, 0, 0))
    img.save(path)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["full", "tools_only"], required=True)
    a = ap.parse_args()
    proj = ROOT / "arena" / f"_probe_wb_{a.mode}"
    proj.mkdir(parents=True, exist_ok=True)
    (proj / "inputs").mkdir(exist_ok=True)
    make_picture(proj / "inputs" / "probe.png")
    out = ROOT / "control" / "environment" / f"probe_wb_{a.mode}"
    out.mkdir(parents=True, exist_ok=True)
    sink = Sink(out / "events.jsonl")
    wb = Workbench()
    rep: dict = {"mode": a.mode, "project": str(proj), "started": time.strftime("%Y-%m-%dT%H:%M:%S")}

    t0 = time.monotonic()
    rep["open"] = wb.open_project(str(proj))
    rep["settings_applied"] = wb.set_settings(str(proj), permission_mode="auto",
                                              settings={**CONDITION_SETTINGS, "knowledge_mode": a.mode})
    rep["settings_readback"] = wb.get_settings(str(proj))
    rep["open_seconds"] = round(time.monotonic() - t0, 1)

    # warm-up turn (the same one every trial gets; not scientific)
    t1 = time.monotonic()
    s = wb.send(str(proj), WARMUP, title="warmup")
    thread = s["thread_id"]
    st = wb.follow(thread, str(proj), sink, deadline=time.monotonic() + 900)
    rep["warmup"] = {"thread": thread, "ended": st["ended"], "seconds": round(time.monotonic() - t1, 1), "cursor": st["cursor"]}

    # MCP readiness: the tool list is the frozen schema of this condition
    t2 = time.monotonic()
    mcp = wb.mcp_status(str(proj))
    rep["mcp_seconds"] = round(time.monotonic() - t2, 1)
    (out / "mcp_status.json").write_text(json.dumps(mcp, ensure_ascii=False, indent=1), encoding="utf-8")
    tools = mcp.get("tools") or mcp.get("tool_names") or []
    names = [t.get("name") if isinstance(t, dict) else t for t in tools]
    rep["mcp"] = {"keys": sorted(mcp.keys()), "n_tools": len(names), "tools": sorted(n for n in names if n)}

    # image check: does the attachment reach the model?
    t3 = time.monotonic()
    s2 = wb.send(str(proj), IMAGE_Q, thread_id=thread,
                 attachments=[{"rel": "inputs/probe.png", "kind": "image", "name": "probe.png"}])
    st2 = wb.follow(thread, str(proj), sink, cursor=st["cursor"], deadline=time.monotonic() + 600)
    rep["image_turn"] = {"ended": st2["ended"], "seconds": round(time.monotonic() - t3, 1), "cursor": st2["cursor"]}
    # pull the assistant text of that turn from the transcript
    tr = wb.transcript(thread, str(proj), limit=200)
    texts = []
    for ev in tr.get("events", []):
        if ev.get("kind") in ("agent_message", "assistant_message", "message") or "text" in ev:
            texts.append(str(ev.get("text") or ev.get("content") or "")[:300])
    rep["image_turn"]["assistant_texts_tail"] = texts[-3:]
    rep["transcript_kinds"] = sorted({e.get("kind") for e in tr.get("events", []) if isinstance(e, dict)})
    rep["agents_md_marker"] = (proj / "AGENTS.md").read_text(encoding="utf-8").splitlines()[0] if (proj / "AGENTS.md").exists() else None
    rep["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    (out / "probe_report.json").write_text(json.dumps(rep, ensure_ascii=False, indent=2), encoding="utf-8")
    sink.close()
    print(json.dumps({k: rep[k] for k in ("open_seconds", "warmup", "mcp_seconds", "image_turn", "agents_md_marker")}, ensure_ascii=False))
    print("tools:", rep["mcp"]["n_tools"], rep["mcp"]["tools"][:60])
    return 0


if __name__ == "__main__":
    sys.exit(main())
