"""Match the primitive Dirichlet theorem and its independent verification records by digest."""
import argparse
from hashlib import sha256
from pathlib import Path
import subprocess

SOURCE = "https://raw.githubusercontent.com/leanprover-community/riemann-hypothesis/371a1a26d5e17f604428e3a8ce6ce341cf517434/src/"
RECORDS = "https://math.mit.edu/~drew/grh/"
FILES = (
    ("Main.lean", SOURCE, "017645ca448f2fd1f6fd5c4c2cb2c6e96e472725be13f6e180ca6e1574360939"),
    ("Statements.lean", SOURCE, "79ef7ca776e42abe1f46caee13ea09bd19f1e2d0b1620d8d5d200258f33b47bf"),
    ("audit.log", RECORDS, "0a1473fce6eed840011f4a76b1c51558fc6f03edd1e6bd1ee542e6c738d87a62"),
    ("con-leche-compare.log", RECORDS, "5a068bf55f24904eea9f05eda947c63d43bdee3c6ab0490e43943a3bd9af1f6e"),
    ("nanoda.log", RECORDS, "d4c5360b02d8dd885150c1cc28b616cb6fb5bf5b20d142558d79a7b15ae9d9df"),
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", type=Path,
                        help="optional local directory containing the five files with published names")
    path = parser.parse_args().directory
    files = {}
    for name, base, expected in FILES:
        data = ((path/name).read_bytes() if path else subprocess.check_output(
            ["curl", "-fsSL", "--retry", "3", "--max-time", "60", base+name]))
        if sha256(data).hexdigest() != expected:
            raise ValueError("Digest mismatch for " + name)
        files[name] = data.decode()
    theorem = "LFunction_eq_zero_re_of_isPrimitive"
    if not all((theorem + " : PrimitiveDirichletLZerosOnLine := critAll" in files["Main.lean"],
                "χ.IsPrimitive →" in files["Statements.lean"],
                "0 < s.re → s.re < 1 → DirichletCharacter.LFunction χ s = 0 → s.re = 1 / 2" in files["Statements.lean"],
                f"'{theorem}' depends on axioms: [propext, Classical.choice, Quot.sound]" in files["audit.log"],
                theorem+" justified by the solution" in files["con-leche-compare.log"],
                "Exit status: 0" in files["nanoda.log"])):
        raise ValueError("Unexpected upstream statement or verification record")
    print("5 pinned upstream files match: primitive Dirichlet strip theorem, three standard axioms, two external kernel records")


if __name__ == "__main__":
    main()
