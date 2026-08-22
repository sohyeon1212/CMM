from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

from cmm.cli import build_parser, main


def test_cli_exposes_production_and_report_commands() -> None:
    help_text = build_parser().format_help()
    assert "production-targets" in help_text
    assert "report" in help_text


def test_report_validate_returns_nonzero_and_machine_readable_errors(
    tmp_path, capsys
) -> None:
    missing = tmp_path / "not-a-run"
    assert main(["report", "validate", str(missing), "--json"]) == 1
    payload = json.loads(capsys.readouterr().out)
    assert payload["valid"] is False
    assert payload["issues"]


def test_production_targets_analysis_only_uses_configured_run_directory(
    tmp_path, monkeypatch, capsys
) -> None:
    config_path = tmp_path / "config.json"
    output = tmp_path / "run"
    config_path.write_text(
        json.dumps(
            {
                "model_path": "model.xml",
                "product": "EX_product_e",
                "output_dir": str(output),
                "run_single_knockout": False,
                "run_strain_design": False,
                "run_amplification": False,
                "validation": {"enabled": False},
            }
        ),
        encoding="utf-8",
    )

    import cmm.workflows

    seen: list[Path] = []

    def fake_run(config):
        seen.append(config.output_dir)
        return SimpleNamespace(run_directory=output.resolve())

    monkeypatch.setattr(cmm.workflows, "run_production_target_discovery", fake_run)
    assert (
        main(
            [
                "production-targets",
                "--config",
                str(config_path),
                "--analysis-only",
            ]
        )
        == 0
    )
    assert seen == [output]
    assert capsys.readouterr().out.strip() == str(output.resolve())


def test_production_targets_requires_an_output_directory(tmp_path, capsys) -> None:
    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps({"model_path": "model.xml", "product": "EX_product_e"}),
        encoding="utf-8",
    )
    assert main(["production-targets", "--config", str(config_path)]) == 1
    assert "config.output_dir" in capsys.readouterr().err
