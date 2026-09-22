"""Remove the two test cases' fingerprints from the frozen knowledge copy used by condition H1.

Only lines that name NU-1000 / the 'hex' development case or carry its numeric fingerprint
(cell volume, R1 values of earlier runs) are rewritten; the general rule each line teaches is kept.
Writes provenance/knowledge_redactions.json with before/after SHA-256 of every touched file.
"""
import hashlib, json, sys
from pathlib import Path

ROOT = Path(r"H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921")
KN = ROOT / "engine" / "CrystalPilot" / "knowledge"
OUT = ROOT / "control" / "provenance" / "knowledge_redactions.json"

EDITS = [
    ("skills/element-assignment-audit.md",
     "是最后一步而不是依据。pa1 hex-l1-r1 把 NU-1000 的 Zr 全部标成 Zn\n  （R1 0.0817 vs 0.0904），validate_structure 的 \"Zn CN=8 outside [4,6]\"",
     "是最后一步而不是依据。一次批测把某 Zr 簇框架的 Zr 全部标成 Zn\n  （改判后 R1 反而更低），validate_structure 的 \"Zn CN=8 outside [4,6]\"",
     "names the NU-1000 case, its element identity and the R1 values of earlier runs"),
    ("skills/framework-solve-ladder.md",
     "  pa1 实测：hex（22 000 Å³、6/mmm）定相 150–840 s、群搜索 300 s、\n  元素指认 90 s；cage（19 000 Å³、2/m）定相 110–830 s。",
     "  pa1 实测：大胞高对称框架定相 150–840 s、群搜索 300 s、\n  元素指认 90 s；cage（19 000 Å³、2/m）定相 110–830 s。",
     "cell volume and Laue class fingerprint the NU-1000 development case"),
    ("skills/data-ingest-space-group-protocol.md",
     "逐个调用，更不要用 change_space_group 逐群声明来\"看 Rint\"（pa1 hex\n  两个 run 各花 13 min 在这上面）。",
     "逐个调用，更不要用 change_space_group 逐群声明来\"看 Rint\"（有的\n  run 各花 13 min 在这上面）。",
     "names the 'hex' development case (NU-1000)"),
]

def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()

records = []
for rel, old, new, why in EDITS:
    p = KN / rel
    text = p.read_text(encoding="utf-8")
    if old not in text:
        print("NOT FOUND:", rel); sys.exit(1)
    before = sha(p)
    p.write_text(text.replace(old, new, 1), encoding="utf-8")
    records.append({"file": rel, "sha256_before": before, "sha256_after": sha(p), "reason": why,
                    "removed_text": old, "replacement_text": new})
    print("redacted", rel)
# a second pass proves no case fingerprint survives anywhere in the snapshot
needles = ["NU-1000", "nu1000", "NU1000", "22 000 Å", "hex-l1", "pa1 hex", "alanine", "丙氨酸", "2360269", "pgw240033"]
left = []
for f in KN.rglob("*.md"):
    t = f.read_text(encoding="utf-8")
    for n in needles:
        if n.lower() in t.lower():
            left.append({"file": str(f.relative_to(KN)), "needle": n})
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps({"edits": records, "residual_matches": left, "needles": needles}, ensure_ascii=False, indent=2), encoding="utf-8")
print("residual matches:", left)
