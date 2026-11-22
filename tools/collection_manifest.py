"""Write or check the checksum list for tracked files outside the paper manifests."""

from pathlib import Path
import argparse
import hashlib
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    manifest = root / "COLLECTION.sha256"
    tracked = subprocess.check_output(["git", "ls-files", "-z"], cwd=root)
    paths = {Path(name.decode()) for name in tracked.split(b"\0") if name}
    covered = {Path(manifest.name)}
    for paper in sorted(path for path in paths if path.name == "SHA256SUMS"):
        covered.add(paper)
        for line in (root / paper).read_text().splitlines():
            _, filename = line.split("  ", 1)
            covered.add(paper.parent / filename)
    listed = sorted(paths - covered, key=lambda path: path.as_posix())
    if args.write:
        entries = []
        for path in listed:
            digest = hashlib.sha256((root / path).read_bytes()).hexdigest()
            entries.append(f"{digest}  {path.as_posix()}\n")
        manifest.write_text("".join(entries))
    else:
        actual = [Path(line.split("  ", 1)[1]) for line in manifest.read_text().splitlines()]
        if actual != listed:
            raise SystemExit("Collection file list differs from the tracked files. Check the paper manifests and run python3 tools/collection_manifest.py --write.")


if __name__ == "__main__":
    main()
