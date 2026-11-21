"""Compile the variable-width kernel lemmas, check their axioms and replay their kernel proofs."""
import os
from pathlib import Path
import subprocess
import tempfile

COMPILER = "819816b2e0a3bf405af45ae5c7af2491d8f5bee6"
THEOREMS = ("triangle", "hyperbolic_stencil", "extremal_weight", "subtract_background", "single_power")
ALLOWED = {"propext", "Classical.choice", "Quot.sound"}


def main():
    repo = Path(__file__).resolve().parent
    if COMPILER not in subprocess.check_output(["lean", "--version"], text=True):
        raise RuntimeError("Use the compiler pinned in lean-toolchain.")
    with tempfile.TemporaryDirectory() as directory:
        obj = Path(directory) / "WidthCusp.olean"
        compiler = subprocess.run(
            ["lean", "-R", str(repo), "-o", str(obj), str(repo / "WidthCusp.lean")],
            text=True, stdout=subprocess.PIPE, check=True)
        lines = compiler.stdout.strip().splitlines()
        if len(lines) != len(THEOREMS):
            raise RuntimeError(compiler.stdout)
        for theorem, line in zip(THEOREMS, lines):
            prefix = f"'WidthCusp.{theorem}' depends on axioms: ["
            if not line.startswith(prefix) or not line.endswith("]"):
                raise RuntimeError(line)
            if set(line[len(prefix):-1].split(", ")) != ALLOWED:
                raise RuntimeError(line)
        env = os.environ.copy()
        env["LEAN_PATH"] = directory + os.pathsep + env.get("LEAN_PATH", "")
        subprocess.run(["leanchecker", "WidthCusp"], env=env, check=True,
                       stdout=subprocess.DEVNULL)
    print(f"Lean 4.33.1 at {COMPILER}; kernel replay passed for WidthCusp")
    for line in lines:
        print(line)


if __name__ == "__main__":
    main()
