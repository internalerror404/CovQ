"""J8: SHA-256 manifest over source and results, plus an environment lock."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TRACKED = ("prototype/src", "prototype/tests", "prototype/scripts", "results",
           "docs", "charter", "STATUS.md", "REGISTRATION.md", "README.md",
           "prototype/README.md", "prototype/requirements.txt")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git(*args, default="unknown"):
    try:
        return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True,
                              text=True, check=True).stdout.strip()
    except Exception:
        return default


def main() -> int:
    files = []
    for rel in TRACKED:
        base = ROOT / rel
        if base.is_file():
            files.append(base)
        elif base.is_dir():
            files.extend(p for p in sorted(base.rglob("*"))
                         if p.is_file() and "__pycache__" not in p.parts)
    entries = [{"path": str(p.relative_to(ROOT)), "bytes": p.stat().st_size,
                "sha256": sha256(p)} for p in sorted(files)]

    import numpy
    import scipy
    manifest = {
        "schema_version": "covq.release_manifest/0.4",
        "release_commit": git("rev-parse", "HEAD"),
        # Scoped to source, for the same reason the records are: this script
        # writes into results/, so a whole-tree check would always report dirty
        # and the flag would carry no information.
        "source_dirty": bool(git("status", "--porcelain", "--",
                                 "prototype/src", "prototype/tests",
                                 "prototype/scripts", "paper", "docs", "charter",
                                 "STATUS.md", "REGISTRATION.md")),
        "environment": {
            "python": sys.version.split()[0],
            "numpy": numpy.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
        },
        "n_files": len(entries),
        "files": entries,
    }
    # Record the tag here rather than patching it in afterwards: this generator
    # rebuilds the manifest from scratch, so anything added to the JSON by hand
    # is silently dropped on the next run.
    tag = git("describe", "--tags", "--exact-match", default="")
    manifest["release_tag"] = {
        "name": tag or None,
        "target_commit": git("rev-list", "-n", "1", tag) if tag else None,
        "pushed": False,
        "reason": "this environment's git proxy refuses refs/tags; branch refs "
                  "push normally, so the commit SHA is the durable identifier",
    }
    manifest["manifest_sha256"] = hashlib.sha256(
        json.dumps(entries, sort_keys=True).encode()).hexdigest()
    out = ROOT / "results/RELEASE_MANIFEST.json"
    out.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"{len(entries)} files  manifest_sha256={manifest['manifest_sha256'][:16]}...")
    print(f"release_commit={manifest['release_commit'][:12]}  source_dirty={manifest['source_dirty']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
