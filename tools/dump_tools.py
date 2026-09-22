"""Freeze the MCP tool schemas each workbench condition exposes (names, descriptions, input schemas),
built from the engine copy's registry exactly as the MCP server builds it for a project.

    dump_tools.py --mode full|tools_only --project DIR
"""
import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

ROOT = Path(r"H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["full", "tools_only"], required=True)
    ap.add_argument("--project", required=True)
    a = ap.parse_args()
    if a.mode == "tools_only":
        os.environ["CRYSTALPILOT_KNOWLEDGE_MODE"] = "tools_only"
    sys.path.insert(0, str(ROOT / "engine" / "CrystalPilot"))
    from crystalpilot.refine.registry import refinement_registry  # noqa: E402
    from crystalpilot.refine.project import RefineProject  # noqa: E402
    reg = refinement_registry(RefineProject(Path(a.project)))
    specs = reg.specs()
    out = ROOT / "control" / "freeze" / f"tool_schemas_{a.mode}.json"
    payload = {"mode": a.mode, "n_tools": len(specs), "names": sorted(s.get("name") for s in specs), "specs": specs}
    out.write_text(json.dumps(payload, ensure_ascii=False, indent=1, default=str), encoding="utf-8")
    digest = hashlib.sha256(json.dumps(specs, sort_keys=True, default=str).encode("utf-8")).hexdigest()
    print(a.mode, "tools:", len(specs), "sha256:", digest[:16])
    return 0


if __name__ == "__main__":
    sys.exit(main())
