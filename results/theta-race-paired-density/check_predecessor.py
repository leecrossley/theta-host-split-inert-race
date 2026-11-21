"""Check the exact predecessor transcript and the ordered table of certified zero brackets."""
from hashlib import sha256
from pathlib import Path


def main():
    evidence = Path(__file__).resolve().parent / "evidence"
    files = {
        "paired-window-check.txt": "c454eaa9ebff709db4df3668def47cc4bcba2f4879b7118a2dacc1ffbc093212",
        "zero-brackets.tsv": "e5483b464160fbce34313daa263f18dae8cc22035044a2053784a274e13a10f3",
    }
    for name, digest in files.items():
        if sha256((evidence / name).read_bytes()).hexdigest() != digest:
            raise ValueError("Wrong predecessor file: " + name)
    print("Exact 400H v1.0.0 bracket table and signed-window verification transcript matched.")


if __name__ == "__main__":
    main()
