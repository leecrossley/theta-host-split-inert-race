"""Compile the algebraic-height lemmas, check their axioms and replay their kernel proofs."""
import os
from pathlib import Path
import subprocess
import tempfile

COMPILER = "819816b2e0a3bf405af45ae5c7af2491d8f5bee6"
THEOREMS = ("protected_log", "affine_zeros", "square_refinement", "threshold_cell", "eligible", "simultaneous_cells",
            "inverse_discriminant", "common_normalizer", "repeated_cell", "paired_roots",
            "cell_energy", "energy_jump", "simple_root_obstruction",
            "ratio_relation", "recover_scaled_slope", "fresh_prime_order",
            "dilation_elimination", "dilation_amplitude", "passage_length", "farey_denominator",
            "log_midpoint", "isolate_character", "lattice_margin",
            "signed_vandermonde", "reciprocal_secant", "four_node_quadratic", "two_weight_interval",
            "normalized_first_hit", "first_branch_square", "one_branch_degree", "pair_unit_quartic", "scaled_elliptic_cubic",
            "harmonic_numerator", "half_point_tangent", "half_map_degree", "genus_two_coefficients")
ALLOWED = {"propext", "Classical.choice", "Quot.sound"}


def main():
    repo = Path(__file__).resolve().parent
    if COMPILER not in subprocess.check_output(["lean", "--version"], text=True):
        raise RuntimeError("Use the compiler pinned in lean-toolchain.")
    with tempfile.TemporaryDirectory() as directory:
        obj = Path(directory) / "AlgebraicHeights.olean"
        compiler = subprocess.run(
            ["lean", "-R", str(repo), "-o", str(obj), str(repo / "AlgebraicHeights.lean")],
            text=True, stdout=subprocess.PIPE, check=True)
        lines = compiler.stdout.strip().splitlines()
        if len(lines) != len(THEOREMS):
            raise RuntimeError(compiler.stdout)
        for theorem, line in zip(THEOREMS, lines):
            prefix = f"'AlgebraicHeights.{theorem}' depends on axioms: ["
            if not line.startswith(prefix) or not line.endswith("]"):
                raise RuntimeError(line)
            expected = {"propext", "Quot.sound"} if theorem in {"one_branch_degree", "half_map_degree", "genus_two_coefficients"} else ALLOWED
            if set(line[len(prefix):-1].split(", ")) != expected:
                raise RuntimeError(line)
        env = os.environ.copy()
        env["LEAN_PATH"] = directory + os.pathsep + env.get("LEAN_PATH", "")
        subprocess.run(["leanchecker", "AlgebraicHeights"], env=env, check=True,
                       stdout=subprocess.DEVNULL)
    print(f"Lean 4.33.1 at {COMPILER}; kernel replay passed for AlgebraicHeights")
    for line in lines:
        print(line)


if __name__ == "__main__":
    main()
