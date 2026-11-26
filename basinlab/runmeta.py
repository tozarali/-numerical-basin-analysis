"""Uniform provenance block for every experiment summary.json."""
import json
import platform
import subprocess
import sys
import time
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import scipy


def git_commit(root):
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root,
                             capture_output=True, text=True, timeout=5)
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=root,
                               capture_output=True, text=True, timeout=5)
        if out.returncode:
            return None
        return out.stdout.strip() + ("+dirty" if dirty.stdout.strip() else "")
    except (OSError, subprocess.SubprocessError):
        return None


def environment(root):
    return {"python": sys.version.split()[0], "numpy": np.__version__,
            "scipy": scipy.__version__, "platform": platform.platform(),
            "git_commit": git_commit(root)}


@contextmanager
def experiment(name, config, *, status="exploratory"):
    """Usage:
        with experiment("dsep_exploratory", cfg) as run:
            ...; run["results"] = {...}
    Writes data/<name>/summary.json even if the body raises (status=failed).
    status must be 'exploratory' or 'confirmatory' (confirmatory requires a
    committed protocol/PREREG.md)."""
    if status not in {"exploratory", "confirmatory"}:
        raise ValueError("status must be exploratory or confirmatory")
    root = Path(__file__).resolve().parents[1]
    out = root / "data" / name
    out.mkdir(parents=True, exist_ok=True)
    run = {"name": name, "status": status, "configuration": config,
           "environment": environment(root), "results": {}, "out_dir": str(out)}
    started = time.perf_counter()
    try:
        yield run
        run["outcome"] = "completed"
    except BaseException as error:
        run["outcome"] = f"failed: {type(error).__name__}: {error}"
        raise
    finally:
        run["runtime_seconds"] = round(time.perf_counter() - started, 3)
        record = dict(run, out_dir=out.relative_to(root).as_posix())   # no absolute paths on disk
        (out / "summary.json").write_text(json.dumps(record, indent=2, default=float))
