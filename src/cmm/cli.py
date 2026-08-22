"""Command line entry point for CMM."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
import json
from pathlib import Path
import sys

from cmm import __version__


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="cmm")
    parser.add_argument("--version", action="store_true", help="Print the CMM version.")
    commands = parser.add_subparsers(dest="command")

    production = commands.add_parser(
        "production-targets",
        help="Run the canonical production-target-discovery workflow from JSON config.",
    )
    production.add_argument(
        "--config",
        required=True,
        type=Path,
        help="UTF-8 JSON ProductionWorkflowConfig file.",
    )
    production.add_argument(
        "--analysis-only",
        action="store_true",
        help="Write scientific artifacts without invoking the R publication renderer.",
    )
    production.add_argument(
        "--renderer",
        default="nature-r",
        choices=("nature-r",),
        help="Publication renderer used after a successful analysis.",
    )

    report = commands.add_parser("report", help="Render or validate a schema-v2 run.")
    report_commands = report.add_subparsers(dest="report_command", required=True)
    render = report_commands.add_parser(
        "render", help="Render R figures and linked/standalone HTML reports."
    )
    render.add_argument("run_dir", type=Path)
    render.add_argument("--renderer", default="nature-r", choices=("nature-r",))
    validate = report_commands.add_parser(
        "validate", help="Validate the run manifest and publication source artifacts."
    )
    validate.add_argument("run_dir", type=Path)
    validate.add_argument(
        "--json",
        action="store_true",
        dest="as_json",
        help="Print the validation result as JSON.",
    )
    return parser


def _validation_payload(report) -> dict[str, object]:
    return {
        "valid": report.valid,
        "issues": list(report.issues),
        "warnings": list(report.warnings),
        "run_directory": str(report.run.root) if report.run is not None else None,
    }


def _run_production(args: argparse.Namespace) -> int:
    from cmm.reporting import render_production_report, validate_production_run
    from cmm.workflows import ProductionWorkflowConfig, run_production_target_discovery

    config = ProductionWorkflowConfig.from_json(args.config)
    if config.output_dir is None:
        raise ValueError(
            "production-targets requires config.output_dir so the analysis has a "
            "self-contained run directory"
        )
    result = run_production_target_discovery(config)
    if result.run_directory is None:  # guarded above, retained as an invariant check
        raise RuntimeError(
            "production workflow completed without exporting a run directory"
        )
    if args.analysis_only:
        print(result.run_directory)
        return 0

    bundle = render_production_report(result.run_directory, renderer=args.renderer)
    validation = validate_production_run(result.run_directory)
    validation.raise_for_errors()
    print(
        json.dumps(
            {
                "run_directory": str(result.run_directory),
                "report_html": str(bundle.report.report_html),
                "report_standalone_html": str(bundle.report.report_standalone_html),
                "figure_manifest": str(bundle.figures.path),
                "valid": validation.valid,
            },
            indent=2,
        )
    )
    return 0


def _run_report(args: argparse.Namespace) -> int:
    from cmm.reporting import render_production_report, validate_production_run

    if args.report_command == "render":
        bundle = render_production_report(args.run_dir, renderer=args.renderer)
        print(bundle.report.report_html)
        print(bundle.report.report_standalone_html)
        return 0

    validation = validate_production_run(args.run_dir)
    payload = _validation_payload(validation)
    if args.as_json:
        print(json.dumps(payload, indent=2))
    elif validation.valid:
        validated = validation.raise_for_errors()
        print(f"valid schema-v2 production run: {validated.root}")
        for warning in validation.warnings:
            print(f"warning: {warning}")
    else:
        for issue in validation.issues:
            print(f"error: {issue}", file=sys.stderr)
    return 0 if validation.valid else 1


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.version:
        print(__version__)
        return 0

    try:
        if args.command == "production-targets":
            return _run_production(args)
        if args.command == "report":
            return _run_report(args)
    except Exception as error:
        print(f"cmm: {type(error).__name__}: {error}", file=sys.stderr)
        return 1

    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
