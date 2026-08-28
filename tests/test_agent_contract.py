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


TRANSFORMATION_SKILL_ROOT = (
    ROOT / ".agents" / "skills" / "cmm-transformation-engineering"
)


def test_transformation_skill_metadata_and_public_boundaries_are_current() -> None:
    from cmm.workflows.transformation import (
        TransformationWorkflowConfig,
        TransformationWorkflowResult,
        run_transformation_target_discovery,
    )

    skill = (TRANSFORMATION_SKILL_ROOT / "SKILL.md").read_text(encoding="utf-8")
    assert skill.startswith("---\n")
    frontmatter = yaml.safe_load(skill.split("---", 2)[1])
    assert frontmatter["name"] == "cmm-transformation-engineering"
    # The description is how a host decides not to route a production request here.
    assert "cmm-production-engineering" in frontmatter["description"]

    openai = yaml.safe_load(
        (TRANSFORMATION_SKILL_ROOT / "agents" / "openai.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert "$cmm-transformation-engineering" in openai["interface"]["default_prompt"]
    assert openai["policy"]["allow_implicit_invocation"] is True

    parser = build_parser()
    assert "transformation-targets" in parser.format_help()
    parsed = parser.parse_args(["transformation-targets", "--config", "run.json"])
    assert parsed.config == Path("run.json")

    assert TransformationWorkflowConfig is not None
    assert TransformationWorkflowResult is not None
    assert callable(run_transformation_target_discovery)


def test_transformation_skill_relative_markdown_links_resolve() -> None:
    markdown_files = [
        TRANSFORMATION_SKILL_ROOT / "SKILL.md",
        *TRANSFORMATION_SKILL_ROOT.glob("references/*.md"),
    ]
    for markdown in markdown_files:
        text = markdown.read_text(encoding="utf-8")
        for destination in re.findall(r"\[[^]]+\]\(([^)]+)\)", text):
            if "://" in destination or destination.startswith("#"):
                continue
            path = destination.split("#", 1)[0]
            assert (markdown.parent / path).resolve().exists(), (
                f"broken relative link in {markdown.relative_to(ROOT)}: {destination}"
            )


def test_transformation_skill_states_the_deviations_it_must_disclose() -> None:
    # These three are structural, not incidental: a run that omits them reads as reproducing a
    # published pipeline it does not reproduce.
    skill = (TRANSFORMATION_SKILL_ROOT / "SKILL.md").read_text(encoding="utf-8")
    assert "iMAT" in skill
    assert "chosen, not derived" in skill  # epsilon
    assert "full, not partial" in skill  # coupling
    # And the one input no inspection can check.
    assert "Never infer it from file names" in skill


def test_agents_md_names_both_skills_and_their_scenarios_exist() -> None:
    # There is no scenario router any more: the skill is the routing surface, and AGENTS.md
    # only has to name the two skills so a host that reads it up front lands in the same place
    # as one that matched a skill description.
    agents = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    for skill_name in ("cmm-production-engineering", "cmm-transformation-engineering"):
        assert skill_name in agents, f"{skill_name} is missing from AGENTS.md"
    assert "the skill is the execution contract" in agents

    scenarios = ROOT / "docs" / "scenarios"
    assert (scenarios / "SC-01-production-target-discovery.md").is_file()
    assert (scenarios / "SC-04-transformation-target-discovery.md").is_file()
    frontmatter = yaml.safe_load(
        (scenarios / "SC-04-transformation-target-discovery.md")
        .read_text(encoding="utf-8")
        .split("---", 2)[1]
    )
    assert frontmatter["id"] == "SC-04"
    assert "source_expression" in frontmatter["requires"]
    assert "target_expression" in frontmatter["requires"]


def test_each_installed_workflow_is_runnable_from_its_skill_alone() -> None:
    """The skill is the execution contract; the scenario is reference material.

    A host that fires a skill loads only SKILL.md. If the command, the entry point, or a
    blocking capability gate lived in the scenario document instead, the run would either fail
    or silently skip a check — so each skill has to carry them itself.
    """

    for root, command, entry_point, gate in (
        (
            SKILL_ROOT,
            "cmm production-targets",
            "run_production_target_discovery",
            "MILP",
        ),
        (
            TRANSFORMATION_SKILL_ROOT,
            "cmm transformation-targets",
            "run_transformation_target_discovery",
            "MIQP",
        ),
    ):
        skill = (root / "SKILL.md").read_text(encoding="utf-8")
        assert command in skill, f"{root.name} does not name its CLI command"
        assert entry_point in skill, f"{root.name} does not name its Python entry point"
        assert gate in skill, f"{root.name} does not state its solver gate"
        # Scenario documents are linked as references, never as a prerequisite step.
        assert "docs/scenarios/" in skill or "SC-0" in skill


def test_router_sends_installed_workflows_to_their_skill_not_their_scenario() -> None:
    # Both hosts must converge on the same path: a Claude Code session reads AGENTS.md up
    # front, another host matches the skill description. The router has to name the skill for
    # the two to agree.
    agents = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    for skill_name in ("cmm-production-engineering", "cmm-transformation-engineering"):
        assert skill_name in agents, f"{skill_name} is missing from the router"
    assert "the skill is the execution contract" in agents
