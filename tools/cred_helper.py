"""apiKeyHelper for the isolated Claude Code home used by condition P1.

Prints the gateway bearer token read from the product's credential file, so the key never sits in a
settings file, an environment variable or a transcript. Same parsing rule as the product's own hook.
"""
import re
import sys
from pathlib import Path

CRED = Path(r"H:\CrystalPilot\secrets\crystalpilot.txt")


def main() -> int:
    text = CRED.read_text(encoding="utf-8")
    m = re.search(r"Key[：:]\s*(\S+)", text)
    key = m.group(1) if m else text.strip()
    if not key or any(ch.isspace() for ch in key):
        return 1
    sys.stdout.write(key)
    return 0


if __name__ == "__main__":
    sys.exit(main())
