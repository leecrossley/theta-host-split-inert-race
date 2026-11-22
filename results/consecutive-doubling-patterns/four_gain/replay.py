"""Replay the signed octave counts and the higher-degree prime-ideal contributions."""

import argparse
import gzip
import hashlib
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

from check_octaves import EVIDENCE, OCTAVES

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE / "primesieve-v12.18.tar.gz"
SOURCE_SHA256 = "3bb3afdfb891f9e9bcc9e10ee24de5fa2e886277f30570c43f5350ecfd45b1a9"


def build_programs(work, threads, rebuild):
    if hashlib.sha256(ARCHIVE.read_bytes()).hexdigest() != SOURCE_SHA256:
        raise ValueError("The primesieve source archive has changed")
    subprocess.run(["tar", "-xzf", ARCHIVE, "-C", work], check=True)
    source = work / "kimwalisch-primesieve-750cf51"
    cxx = os.environ.get("CXX", "c++")
    flags = [cxx, "-std=c++17", "-O3", "-pthread"]
    sieve_flags = [
        *flags,
        "-DNDEBUG",
        "-I" + str(source / "include"),
        "-I" + str(source / "src"),
    ]

    def compile_object(path):
        obj = work / (path.stem + ".o")
        subprocess.run([*sieve_flags, "-c", path, "-o", obj], check=True)
        return obj

    with ThreadPoolExecutor(min(threads, 8)) as pool:
        objects = list(pool.map(compile_object, sorted((source / "src").glob("*.cpp"))))
    subprocess.run(
        [*flags, HERE / "recheck_cells.cpp", "-o", work / "cells"], check=True
    )
    subprocess.run(
        [
            *sieve_flags,
            HERE / "higher_degrees.cpp",
            *objects,
            "-o",
            work / "degrees-primary",
        ],
        check=True,
    )
    subprocess.run(
        [
            *flags,
            "-DINDEPENDENT",
            HERE / "higher_degrees.cpp",
            "-o",
            work / "degrees-independent",
        ],
        check=True,
    )
    if rebuild:
        subprocess.run(
            [
                *sieve_flags,
                HERE / "octave_moments.cpp",
                *objects,
                "-o",
                work / "octaves",
            ],
            check=True,
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rebuild",
        action="store_true",
        help="Recount every prime in all eight stored octaves.",
    )
    parser.add_argument("--threads", type=int, default=min(8, os.cpu_count() or 1))
    args = parser.parse_args()
    if not 1 <= args.threads <= 64:
        parser.error("--threads must be between 1 and 64")

    receipts = EVIDENCE / "four-gain"
    with TemporaryDirectory(prefix="theta-four-gains-") as temp:
        work = Path(temp)
        build_programs(work, args.threads, args.rebuild)
        if args.rebuild:
            for k in OCTAVES:
                computed = subprocess.check_output(
                    [work / "octaves", str(k), str(args.threads)]
                )
                with gzip.open(receipts / f"octave-{k}.txt.gz", "rb") as stream:
                    if computed != stream.read():
                        raise ValueError(
                            f"The regenerated integer moments differ in octave {k}"
                        )
                print(
                    f"Regenerated every one of the 32,768 cells in octave {k}.",
                    file=sys.stderr,
                    flush=True,
                )

        env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}
        transcript = b"".join(
            subprocess.check_output([sys.executable, HERE / name, *paths], env=env)
            for name, paths in (
                ("check_octaves.py", []),
                ("independent_checks.py", [work / "cells"]),
                (
                    "higher_degrees.py",
                    [work / "degrees-primary", work / "degrees-independent"],
                ),
            )
        )
    if transcript != (receipts / "check.txt").read_bytes():
        sys.stderr.buffer.write(transcript)
        raise ValueError("The replay differs from the archived transcript")
    sys.stdout.buffer.write(transcript)


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    main()
