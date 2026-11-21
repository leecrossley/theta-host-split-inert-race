"""Build the transitive imports of ThetaPrimeTrace under the pinned Mathlib environment."""

import argparse
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from graphlib import TopologicalSorter
import os
from pathlib import Path
import re
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="the upstream lean/src directory")
    parser.add_argument("build", type=Path, help="an empty build directory, or a cache from the same pinned sources and toolchain")
    parser.add_argument("--jobs", type=int, default=3)
    parser.set_defaults(target="ThetaPrimeTrace")
    args = parser.parse_args()
    source = args.source.resolve()
    build = args.build.resolve()
    build.mkdir(parents=True, exist_ok=True)
    local = Path(__file__).resolve().parent
    paths = {}
    imports = {}

    def visit(module):
        if module in paths:
            return
        relative = Path(*module.split(".")).with_suffix(".lean")
        path = source / relative
        if not path.is_file():
            path = local / relative
        if not path.is_file():
            return
        paths[module] = path
        imported = re.findall(r"^import\s+([\w.]+)\s*$", path.read_text(), re.MULTILINE)
        for dep in imported:
            visit(dep)
        imports[module] = [dep for dep in imported if dep in paths]

    visit(args.target)
    graph = TopologicalSorter(imports)
    graph.prepare()
    env = os.environ.copy()
    env["LEAN_PATH"] = os.pathsep.join([str(build), *(str(Path(p).resolve())
        for p in env.get("LEAN_PATH", "").split(os.pathsep) if p)])

    def compile_module(module):
        path = paths[module]
        output = build.joinpath(*module.split(".")).with_suffix(".olean")
        requirements = [path, *(build.joinpath(*dep.split(".")).with_suffix(".olean")
                                for dep in imports[module])]
        if output.exists() and all(p.stat().st_mtime <= output.stat().st_mtime for p in requirements):
            return
        output.parent.mkdir(parents=True, exist_ok=True)
        root = source if path.is_relative_to(source) else local
        print("Checking", module, flush=True)
        options = ["--tstack=400000", "-DautoImplicit=false", "-DrelaxedAutoImplicit=false",
                   "-DmaxSynthPendingDepth=3", "-DmaxHeartbeats=4000000", "-Dpp.unicode.fun=true"]
        if module.startswith(("GroundTheory.LocalFields.", "GroundTheory.GlobalFields.")):
            options.append("-DsynthInstance.maxHeartbeats=400000")
        subprocess.run(["lean", *options, "-j2", "-R", str(root), "-o", str(output), str(path)],
                       check=True, cwd=root, env=env)

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        running = {}
        while graph.is_active():
            for module in graph.get_ready():
                running[pool.submit(compile_module, module)] = module
            finished, _ = wait(running, return_when=FIRST_COMPLETED)
            for result in finished:
                result.result()
                graph.done(running.pop(result))
    upstream = sum(path.is_relative_to(source) for path in paths.values())
    print(f"Checked {upstream} upstream and {len(paths) - upstream} local modules for {args.target}.",
          flush=True)


if __name__ == "__main__":
    main()
