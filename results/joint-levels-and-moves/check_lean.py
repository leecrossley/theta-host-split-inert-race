"""Compile the joint level and move lemmas, check their axioms and replay their kernel proofs."""
import os
from pathlib import Path
import subprocess
import tempfile

COMPILER = "819816b2e0a3bf405af45ae5c7af2491d8f5bee6"
THEOREMS = ("projected_norm", "disk_support", "ellipse_support", "projection_gap", "one_sided_pair", "ellipse_axes", "deadline_support", "mixed_stop_band", "contrast_margin", "cover_endpoint",
            "translated_det", "near_perpendicular", "transverse", "simplex_fixed_displacement", "three_point_face", "quartic_elimination", "recover_weight", "positive_mass_sign", "four_sheet_group", "projection_discriminant")
ALLOWED = {"propext", "Classical.choice", "Quot.sound"}


def main():
    repo = Path(__file__).resolve().parent
    if COMPILER not in subprocess.check_output(["lean", "--version"], text=True):
        raise RuntimeError("Use the compiler pinned in lean-toolchain.")
    with tempfile.TemporaryDirectory() as directory:
        obj = Path(directory) / "JointOscillation.olean"
        compiler = subprocess.run(
            ["lean", "-R", str(repo), "-o", str(obj), str(repo / "JointOscillation.lean")],
            text=True, stdout=subprocess.PIPE, check=True)
        lines = compiler.stdout.strip().splitlines()
        if len(lines) != len(THEOREMS):
            raise RuntimeError(compiler.stdout)
        for theorem, line in zip(THEOREMS, lines):
            prefix = f"'JointOscillation.{theorem}' depends on axioms: ["
            if not line.startswith(prefix) or not line.endswith("]"):
                raise RuntimeError(line)
            expected = ALLOWED
            if set(line[len(prefix):-1].split(", ")) != expected:
                raise RuntimeError(line)
        env = os.environ.copy()
        env["LEAN_PATH"] = directory + os.pathsep + env.get("LEAN_PATH", "")
        subprocess.run(["leanchecker", "JointOscillation"], env=env, check=True,
                       stdout=subprocess.DEVNULL)
    print(f"Lean 4.33.1 at {COMPILER}; kernel replay passed for JointOscillation")
    for line in lines:
        print(line)


if __name__ == "__main__":
    main()
