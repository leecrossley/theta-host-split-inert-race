"""Rebuild the host witnesses and fibre identity against their pinned Lean sources."""
import argparse
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from graphlib import TopologicalSorter
from hashlib import sha256
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile


LOCAL = ("HostPrimes", "SecondField", "FiberTrace", "RaceTrace", "ShiftBound")
THEOREMS = ("HostRace.jump_eq", "HostRace.fiber_trace", "SecondHost.normPrime",
            "SecondHost.field_parameters", "SecondHost.field_degrees", "SecondHost.frame_exists",
            "SecondHost.admissible_exists", "SecondHost.concrete_fiber_trace",
            "ShiftBound.width_floor", "ShiftBound.uniform_short_shift")
TRACE = {"DyadicArithmetic": "868f107ed4be631fe5ab70cba0d922ccd81fdc6662c9d2989ea6163f5e2f116d",
         "ThetaPrimeTrace": "f19c9f3dbf5464dae4fb20a383947586fb7e85ec1a527f50d1bd11e684e59eb3"}
COMPILER = "819816b2e0a3bf405af45ae5c7af2491d8f5bee6"
AXIOMS = {"propext", "Classical.choice", "Quot.sound"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("build", type=Path)
    parser.add_argument("trace", type=Path)
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    source, build, trace = (p.resolve() for p in (args.source, args.build, args.trace))
    repo = Path(__file__).resolve().parent
    if sha256((repo / "FiberTrace.lean").read_bytes()).hexdigest() != "6322d03a469e813e8f25db674f3df4ed36f79afe8d1029118fde87d4ad056343":
        raise SystemExit("The fibre theorem differs from the pinned 571/599 release.")
    if not source.is_dir() or build in (repo, source, trace) or args.jobs < 1:
        parser.error("provide a source checkout and a separate build cache")
    if COMPILER not in subprocess.check_output(["lean", "--version"], text=True):
        raise SystemExit("Use the compiler in lean-toolchain.")
    for module, expected in TRACE.items():
        if sha256((trace / f"{module}.lean").read_bytes()).hexdigest() != expected:
            raise SystemExit(f"Trace source mismatch for {module}")
    env = os.environ.copy()
    env["LEAN_PATH"] = os.pathsep.join([str(build), *(str(Path(p).resolve()) for p in
        env.get("LEAN_PATH", "").split(os.pathsep) if p)])
    graph, paths, hashes = {}, {}, {}

    def visit(module):
        if module in graph:
            return
        root = repo if module in LOCAL else (trace if module in TRACE else source)
        path = root / Path(*module.split(".")).with_suffix(".lean")
        if not path.is_file():
            return  # Library or compiler import from the pinned Lake environment.
        text = path.read_bytes()
        paths[module] = (path, root)
        hashes[module] = sha256(text).hexdigest()
        graph[module] = []
        # Decode first: a bytes pattern's \w is ASCII only and would skip modules such as DataParams₂.
        for name in re.findall(r"^import[ \t]+(\S+)[ \t]*$", text.decode("utf-8"), re.M):
            visit(name)
            if name in paths:
                graph[module].append(name)

    visit("RaceTrace")
    visit("ShiftBound")
    upstream = sorted(set(paths).difference(LOCAL, TRACE))
    digest = sha256("".join(f"{m} {hashes[m]}\n" for m in upstream).encode()).hexdigest()
    if digest != "3777adec61ae40822f242c942681f8713c6bcc604ffc020223378a065bedee2a":
        raise SystemExit("Unexpected upstream source digest.")
    build.mkdir(parents=True, exist_ok=True)

    def output(m):
        return build / Path(*m.split(".")).with_suffix(".olean")

    with tempfile.TemporaryDirectory() as temp:
        def compile_module(m):
            path, root = paths[m]
            out = output(m)
            deps = [output(d) for d in graph[m]]
            if m not in LOCAL and out.is_file() and all(
                    p.stat().st_mtime <= out.stat().st_mtime for p in [path, *deps]):
                return
            out.parent.mkdir(parents=True, exist_ok=True)
            flags = ["--tstack=400000", "-DautoImplicit=false", "-DrelaxedAutoImplicit=false",
                     "-DmaxSynthPendingDepth=3", "-DmaxHeartbeats=4000000", "-Dpp.unicode.fun=true"]
            if m.startswith(("GroundTheory.LocalFields.", "GroundTheory.GlobalFields.")):
                flags.append("-DsynthInstance.maxHeartbeats=400000")
            log = Path(temp) / (m + ".log")
            with log.open("w") as handle:
                rc = subprocess.run(["lean", *flags, "-j2", "-R", root, "-o", out, path],
                                    cwd=root, env=env, stdout=handle, stderr=subprocess.STDOUT).returncode
            if rc:
                raise RuntimeError(f"Lean failed on {m}\n{log.read_text()[-8000:]}")

        order = TopologicalSorter(graph)
        order.prepare()
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            running = {}
            completed = 0
            while order.is_active():
                for m in order.get_ready():
                    running[pool.submit(compile_module, m)] = m
                done, _ = wait(running, return_when=FIRST_COMPLETED)
                for task in done:
                    task.result()
                    order.done(running.pop(task))
                    completed += 1
                    if completed % 100 == 0 or completed == len(graph):
                        print(f"Lean build: {completed}/{len(graph)} modules ready", file=sys.stderr, flush=True)
        # leanchecker imports a separate environment for each target in parallel.
        for module in LOCAL:
            print(f"Kernel replay: {module}", file=sys.stderr, flush=True)
            subprocess.run(["leanchecker", module], env=env, cwd=build, check=True,
                           stdout=subprocess.DEVNULL)
        audit = Path(temp) / "Audit.lean"
        audit.write_text("import RaceTrace\nimport ShiftBound\n" + "".join(
            f"#print axioms {t}\n" for t in THEOREMS))
        log = subprocess.check_output(["lean", audit], env=env, cwd=build, text=True)
        audit_results = {}
        for name in THEOREMS:
            if f"'{name}' does not depend on any axioms" in log:
                audit_results[name] = []
                continue
            pat = re.escape(f"'{name}' depends on axioms:") + r"\s*\[([^\]]*)\]"
            matches = re.findall(pat, log)
            if len(matches) != 1 or not set(matches[0].replace(" ", "").split(",")) <= AXIOMS:
                raise RuntimeError(f"Unexpected axiom report for {name}\n{log}")
            audit_results[name] = [item.strip() for item in matches[0].split(",")]
    print(f"Lean 4.33.1 at {COMPILER}")
    print(f"{len(upstream)} upstream source modules, SHA-256 of sorted module names and hashes {digest}")
    print("Pinned trace source hashes match both files")
    print("Kernel replay passed for " + ", ".join(LOCAL))
    for name in THEOREMS:
        print(f"{name}: {', '.join(audit_results[name]) or 'none'}")


if __name__ == "__main__":
    main()
