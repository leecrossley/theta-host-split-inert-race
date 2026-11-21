"""Compile the doubling pattern lemmas, check their axioms and replay their kernel proofs."""
import os
from pathlib import Path
import subprocess
import tempfile

COMPILER = "819816b2e0a3bf405af45ae5c7af2491d8f5bee6"
THEOREMS = ("sign_preserved", "separated_signs", "move_error", "telescoping",
            "excluded_run", "word_count", "reflect_difference", "cubic_majorant")
ALLOWED = {"propext", "Classical.choice", "Quot.sound"}


def main():
    repo = Path(__file__).resolve().parent
    if COMPILER not in subprocess.check_output(["lean", "--version"], text=True):
        raise RuntimeError("Use the compiler pinned in lean-toolchain.")
    with tempfile.TemporaryDirectory() as directory:
        obj = Path(directory) / "DoublingPatterns.olean"
        compiler = subprocess.run(
            ["lean", "-R", str(repo), "-o", str(obj), str(repo / "DoublingPatterns.lean")],
            text=True, stdout=subprocess.PIPE, check=True)
        lines = compiler.stdout.strip().splitlines()
        if len(lines) != len(THEOREMS):
            raise RuntimeError(compiler.stdout)
        for theorem, line in zip(THEOREMS, lines):
            prefix = f"'DoublingPatterns.{theorem}' depends on axioms: ["
            if not line.startswith(prefix) or not line.endswith("]"):
                raise RuntimeError(line)
            expected = ALLOWED
            if not set(line[len(prefix):-1].split(", ")).issubset(expected):
                raise RuntimeError(line)
        env = os.environ.copy()
        env["LEAN_PATH"] = directory + os.pathsep + env.get("LEAN_PATH", "")
        subprocess.run(["leanchecker", "DoublingPatterns"], env=env, check=True,
                       stdout=subprocess.DEVNULL)
    print(f"Lean 4.33.1 at {COMPILER}; kernel replay passed for DoublingPatterns")
    for line in lines:
        print(line)


if __name__ == "__main__":
    main()
