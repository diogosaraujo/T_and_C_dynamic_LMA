"""Where the analysis products go.

The daily and annual effect tables and the drought-year labels all land in one
place: $TC_RESULTS, a sibling of model_run, set by slurm/config.sh to
<TC_ROOT>/result_summary. They are OUTPUTS summarising the runs, so they do not
belong under input_data, and they cannot live in the repo -- era5_daily.csv is
108 MB against GitHub's 100 MB hard limit, and the GCM version covers roughly
15x as many pairs.

WHY THIS MODULE EXISTS. argparse turns "--out era5_daily.csv" into a path
relative to the working directory, which the submit wrappers set to
preprocessing/. That silently wrote a 108 MB table inside the repo, "git add -A"
swept it into a commit, and the push was declined outright. A bare filename now
resolves under $TC_RESULTS instead.

NO SILENT FALLBACK TO THE WORKING DIRECTORY. If a relative name is given and
$TC_RESULTS is not set, that means the submitting shell never sourced
slurm/config.sh, and the correct response is to say so and stop -- writing to
whatever directory the job happens to start in is how this went wrong the first
time. An absolute path is always honoured as given.
"""
from __future__ import annotations

import os
from pathlib import Path


class NoResultsDir(Exception):
    """$TC_RESULTS is needed to place a relative --out and is not set."""


def results_root() -> Path:
    """The configured analysis-output directory, or raise saying how to set it."""
    root = os.environ.get("TC_RESULTS", "").strip()
    if not root:
        raise NoResultsDir(
            "$TC_RESULTS is not set, so a relative --out has nowhere to go.\n"
            "  Run 'source slurm/config.sh' before sbatch, or pass --out as an "
            "absolute path.\n"
            "  Refusing to write to the working directory: that is what put a "
            "108 MB table\n"
            "  inside the repo and got the push declined."
        )
    return Path(root)


def resolve_out(out: Path | str, *, create: bool = True) -> Path:
    """Absolute path for an --out value; relative names go under $TC_RESULTS."""
    p = Path(out)
    if not p.is_absolute():
        p = results_root() / p
    if create:
        p.parent.mkdir(parents=True, exist_ok=True)
    return p


def figures_root() -> Path:
    """$TC_FIGURES -- rendered figures, a sibling of model_run.

    Separate from $TC_RESULTS because a figure is a different product from the
    table behind it: regenerated freely, copied down on its own, and never
    wanted in the repo. Same refusal as results_root -- an unset variable means
    the shell never sourced slurm/config.sh, and writing PNGs into whatever
    directory the job started in is how they end up committed.
    """
    root = os.environ.get("TC_FIGURES", "").strip()
    if not root:
        raise NoResultsDir(
            "$TC_FIGURES is not set, so figures have nowhere to go.\n"
            "  Run 'source slurm/config.sh' before sbatch, or pass an absolute "
            "--out.\n"
            "  Refusing to write into the working directory."
        )
    return Path(root)


def resolve_figure(out: Path | str, *, create: bool = True) -> Path:
    """Absolute path for a figure output; relative names go under $TC_FIGURES."""
    p = Path(out)
    if not p.is_absolute():
        p = figures_root() / p
    if create:
        p.mkdir(parents=True, exist_ok=True)
    return p
