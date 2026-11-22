"""Compile the dyadic prime-support lemmas, check their axioms and replay their kernel proofs."""
import os
from pathlib import Path
import subprocess
import tempfile

COMPILER = "819816b2e0a3bf405af45ae5c7af2491d8f5bee6"
THEOREMS = ("two_zeros", "root_propagates", "prescribe_root", "silent_zero", "invisible_iff",
            "opposite_linear_terms", "dyadic_sum", "phase_quadratic", "half_shift_recovery",
            "right_cell_visible", "affine_second_difference", "one_prime_second_difference",
            "dependent_column_zero")
ALLOWED = {"propext", "Classical.choice", "Quot.sound"}


def main():
    repo = Path(__file__).resolve().parent
    if COMPILER not in subprocess.check_output(["lean", "--version"], text=True):
        raise RuntimeError("Use the compiler pinned in lean-toolchain.")
    with tempfile.TemporaryDirectory() as directory:
        obj = Path(directory) / "DyadicPrimeSupport.olean"
        compiler = subprocess.run(
            ["lean", "-R", str(repo), "-o", str(obj), str(repo / "DyadicPrimeSupport.lean")],
            text=True, stdout=subprocess.PIPE, check=True)
        lines = compiler.stdout.strip().splitlines()
        if len(lines) != len(THEOREMS):
            raise RuntimeError(compiler.stdout)
        for theorem, line in zip(THEOREMS, lines):
            prefix = f"'DyadicPrimeSupport.{theorem}' depends on axioms: ["
            if line == f"'DyadicPrimeSupport.{theorem}' does not depend on any axioms":
                continue
            if not line.startswith(prefix) or not line.endswith("]"):
                raise RuntimeError(line)
            expected = ALLOWED
            if not set(line[len(prefix):-1].split(", ")).issubset(expected):
                raise RuntimeError(line)
        env = os.environ.copy()
        env["LEAN_PATH"] = directory + os.pathsep + env.get("LEAN_PATH", "")
        subprocess.run(["leanchecker", "DyadicPrimeSupport"], env=env, check=True,
                       stdout=subprocess.DEVNULL)
    print(f"Lean 4.33.1 at {COMPILER}; kernel replay passed for DyadicPrimeSupport")
    for line in lines:
        print(line)


if __name__ == "__main__":
    main()
