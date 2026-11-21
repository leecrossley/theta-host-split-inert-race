"""Find inequalities that GitHub can read as HTML inside displayed mathematics."""

import os
from pathlib import Path
import re
import subprocess
import sys


FENCE = re.compile(r"^ {0,3}(`{3,}|~{3,})(.*)$")
LESS_THAN_LETTER = re.compile(r"<[A-Za-z]")


def broken_lines(text):
    fence = None
    math_fence = False
    dollar_block = False
    for number, line in enumerate(text.splitlines(), 1):
        marker = FENCE.match(line)
        if fence:
            if marker and marker[1][0] == fence[0] and len(marker[1]) >= len(fence) and not marker[2].strip():
                fence = None
                continue
            if not math_fence:
                continue
            expression = line
        elif dollar_block:
            expression, closing, _ = line.partition("$$")
            dollar_block = not closing
        elif marker:
            fence = marker[1]
            math_fence = marker[2].strip() == "math"
            continue
        elif re.match(r"^ {0,3}\$\$", line):
            expression, closing, _ = line.lstrip()[2:].partition("$$")
            dollar_block = not closing
        else:
            continue
        if LESS_THAN_LETTER.search(expression):
            yield number


def main():
    paths = subprocess.check_output(["git", "ls-files", "-z", "--", "*.md"]).split(b"\0")
    failed = False
    for name in filter(None, paths):
        path = Path(os.fsdecode(name))
        for number in broken_lines(path.read_text(encoding="utf-8")):
            print(f"{path}:{number}: GitHub can drop part of this inequality. Write \\lt or put a space after <.", file=sys.stderr)
            failed = True
    return int(failed)


if __name__ == "__main__":
    sys.exit(main())
