"""Check replayable copies against their maintained source, or update them with --sync."""
import argparse
import json
from pathlib import Path
import sys

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--sync", action="store_true", help="copy the maintained source to its listed consumers")
args = parser.parse_args()
different = []
for group in json.loads((root / "shared-files.json").read_text()):
    data = (root / group["source"]).read_bytes()
    for relative in group["copies"]:
        destination = root / relative
        if destination.read_bytes() == data:
            continue
        if args.sync:
            destination.write_bytes(data)
            print(f"Updated {relative}")
        else:
            different.append(f"{relative} differs from {group['source']}")
if different:
    sys.exit("\n".join(different))
