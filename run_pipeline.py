"""Run the project pipeline in dependency order from the repository root."""

import argparse
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent


def _run_stage(name: str, script: Path, extra_args: list[str] | None = None) -> None:
    args = extra_args or []
    cmd = [sys.executable, str(script), *args]

    print(f"\n=== {name} ===", flush=True)
    print(f"Running: {' '.join(cmd)}", flush=True)

    subprocess.run(cmd, cwd=REPO_ROOT, check=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run generator, reporting, and analysis in order from the project root."
        )
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Optional seed override passed only to generator/generate.py.",
    )
    args = parser.parse_args(argv)

    generator_args: list[str] = []
    if args.seed is not None:
        generator_args.extend(["--seed", str(args.seed)])

    stages = [
        ("generator", REPO_ROOT / "generator" / "generate.py", generator_args),
        ("reporting", REPO_ROOT / "reporting" / "run_report.py", []),
        ("analysis", REPO_ROOT / "analysis" / "run_analysis.py", []),
    ]

    try:
        for name, script, stage_args in stages:
            _run_stage(name=name, script=script, extra_args=stage_args)
    except subprocess.CalledProcessError as exc:
        script_name = Path(exc.cmd[1]).name if len(exc.cmd) > 1 else "unknown stage"
        print(
            f"\nPipeline failed while running {script_name} "
            f"(exit code {exc.returncode}).",
            file=sys.stderr,
            flush=True,
        )
        return exc.returncode

    print("\nPipeline completed successfully.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
