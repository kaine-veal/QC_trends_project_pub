"""
Site-specific locations, read from environment variables so no local paths are hard-coded.

    QM_RUN_DIR      Sequencer run folders (one folder per run)
    QM_RESULTS_DIR  Pipeline results, organised by year
    QM_OUTPUT_ROOT  Where monthly archives and plots are written

Unset variables fall back to the relative defaults below.
"""
import os
from pathlib import Path


def _path(var: str, default: str) -> Path:
    return Path(os.environ.get(var, default)).expanduser()


RUN_DIR = _path("QM_RUN_DIR", "data/raw")
RESULTS_DIR = _path("QM_RESULTS_DIR", "data/results")
OUTPUT_ROOT = _path("QM_OUTPUT_ROOT", "qc_trends")


def apptainer_binds() -> list:
    """Bind each configured directory into the container at the same path."""
    args = []
    for p in dict.fromkeys(d.resolve() for d in (RUN_DIR, RESULTS_DIR, OUTPUT_ROOT)):
        args += ["--bind", f"{p}:{p}"]
    return args
