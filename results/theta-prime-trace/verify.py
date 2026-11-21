"""Rebuild this repository's modules, replay their declarations, and audit the named theorems."""
import argparse
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

TARGET = 'ThetaPrimeTrace'
MODULES = ['DyadicArithmetic', 'ThetaPrimeTrace']
THEOREMS = ['theta_step_trace', 'theta_prime_trace_extract', 'theta_band_prime_trace']
LEAN_COMMIT = "819816b2e0a3bf405af45ae5c7af2491d8f5bee6"
ALLOWED_AXIOMS = {"propext", "Classical.choice", "Quot.sound"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="pinned upstream lean/src, or an empty directory for a Mathlib-only build")
    parser.add_argument("build", type=Path, help="build directory dedicated to this verifier")
    parser.add_argument("--jobs", type=int, default=3)
    args = parser.parse_args()
    source, build = args.source.resolve(), args.build.resolve()
    if not source.is_dir() or args.jobs < 1:
        parser.error("source must be a directory and jobs must be positive")
    version = subprocess.check_output(["lean", "--version"], text=True).strip()
    if LEAN_COMMIT not in version:
        raise SystemExit("The installed Lean compiler does not match lean-toolchain: " + version)
    repo = Path(__file__).resolve().parent
    build.mkdir(parents=True, exist_ok=True)
    if build == repo or build == source:
        parser.error("use a separate directory for the compiled modules")
    # Force fresh elaboration of this repository's declarations on every invocation.
    for module in MODULES:
        (build / (module + ".olean")).unlink(missing_ok=True)
    subprocess.run([sys.executable, str(repo / "replay.py"), str(source), str(build),
                    "--jobs", str(args.jobs)], check=True)
    env = os.environ.copy()
    paths = [str(Path(p).resolve()) for p in env.get("LEAN_PATH", "").split(os.pathsep) if p]
    env["LEAN_PATH"] = os.pathsep.join([str(build), *paths])
    subprocess.run(["leanchecker", "-v", *MODULES], env=env, check=True)
    with tempfile.TemporaryDirectory() as scratch:
        audit = Path(scratch) / "Audit.lean"
        audit.write_text("import " + TARGET + "\n" + "".join(
            "#print axioms ThetaTransfer." + t + "\n" for t in THEOREMS))
        output = subprocess.check_output(["lean", str(audit)], env=env, text=True)
    for theorem in THEOREMS:
        match = re.search(r"'ThetaTransfer\." + re.escape(theorem) +
                          r"' depends on axioms:\s*\[([^\]]*)\]", output)
        if not match or set(match[1].replace(" ", "").split(",")) != ALLOWED_AXIOMS:
            raise SystemExit("Axiom audit failed for " + theorem + "\n" + output)
    print("Compiler commit:", LEAN_COMMIT)
    print("Kernel replay passed:", ", ".join(MODULES))
    for theorem in THEOREMS:
        print("Axioms for ThetaTransfer." + theorem + ": propext, Classical.choice, Quot.sound")


if __name__ == "__main__":
    main()
