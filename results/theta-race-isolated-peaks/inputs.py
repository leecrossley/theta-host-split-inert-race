"""Read the signed zero brackets and run the bounded integer scan selected by the caller."""
from hashlib import sha256
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parent


def brackets():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    filename = ROOT/"evidence/zero-brackets.tsv"
    assert sha256(filename.read_bytes()).hexdigest() == "e5483b464160fbce34313daa263f18dae8cc22035044a2053784a274e13a10f3"
    transcript = ROOT/"evidence/bracket-check.txt"
    assert sha256(transcript.read_bytes()).hexdigest() == "c454eaa9ebff709db4df3668def47cc4bcba2f4879b7118a2dacc1ffbc093212"
    rows = [tuple(map(int, s.split())) for s in filename.read_text().splitlines()]
    assert len(rows) == 6000 and all(0 < l < u for _, _, l, u in rows)
    ordered = sorted((l, u) for _, _, l, u in rows)
    assert all(u < l for (_, u), (l, _) in zip(ordered, ordered[1:]))
    return rows


def run_scan(source, lines):
    with tempfile.TemporaryDirectory() as directory:
        executable = Path(directory)/"scan"
        subprocess.run([os.environ.get("CC", "cc"), "-std=c11", "-O3", "-Wall", "-Wextra",
                        "-Werror", str(ROOT/source), "-o", str(executable)], check=True)
        subprocess.run([str(executable)], input="\n".join(lines)+"\n", text=True, check=True)
