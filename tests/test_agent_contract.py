from __future__ import annotations

import inspect
from pathlib import Path
import re

import yaml

from cmm.cli import build_parser
from cmm.reporting import render_production_report, validate_production_run
from cmm.workflows import (
    ProductionWorkflowConfig,
    ProductionWorkflowResult,
    run_production_target_discovery,
)


ROOT = Path(__file__).resolve().parents[1]
SKILL_ROOT = ROOT / ".agents" / "skills" / "cmm-production-engineering"


def test_production_skill_metadata_and_public_boundaries_are_current() -> None:
    skill = (SKILL_ROOT / "SKILL.md").read_text(encoding="utf-8")
    assert skill.startswith("---\n")
    frontmatter = yaml.safe_load(skill.split("---", 2)[1])
    assert frontmatter["name"] == "cmm-production-engineering"
    assert "production-target workflow" in frontmatter["description"]

    openai = yaml.safe_load(
        (SKILL_ROOT / "agents" / "openai.yaml").read_text(encoding="utf-8")
    )
    assert "$cmm-production-engineering" in openai["interface"]["default_prompt"]
    assert openai["policy"]["allow_implicit_invocation"] is True
    assert "artifact-validated" in openai["interface"]["default_prompt"]

    parser = build_parser()
    help_text = parser.format_help()
    assert "production-targets" in help_text
    assert "report" in help_text
    production = parser.parse_args(
        ["production-targets", "--config", "workflow.json", "--analysis-only"]
    )
    assert production.analysis_only is True
    assert production.config == Path("workflow.json")
    assert parser.parse_args(["report", "render", "run"]).report_command == "render"
    validation = parser.parse_args(["report", "validate", "run", "--json"])
    assert validation.report_command == "validate"
    assert validation.as_json is True

    assert ProductionWorkflowConfig is not None
    assert ProductionWorkflowResult is not None
    assert callable(run_production_target_discovery)
    assert callable(render_production_report)
    assert callable(validate_production_run)
    renderer = inspect.signature(render_production_report).parameters["renderer"]
    assert renderer.default == "nature-r"


def test_production_skill_relative_markdown_links_resolve() -> None:
    markdown_files = [SKILL_ROOT / "SKILL.md", *SKILL_ROOT.glob("references/*.md")]
    for markdown in markdown_files:
        text = markdown.read_text(encoding="utf-8")
        for destination in re.findall(r"\[[^]]+\]\(([^)]+)\)", text):
            if "://" in destination or destination.startswith("#"):
                continue
            path = destination.split("#", 1)[0]
            assert (markdown.parent / path).resolve().exists(), (
                f"broken relative link in {markdown.relative_to(ROOT)}: {destination}"
            )
