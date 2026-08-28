"""Rendering and validating a transformation run, without needing a solver.

The end-to-end test in ``test_transformation_workflow.py`` needs MIQP and therefore skips
wherever Gurobi or CPLEX is absent — which is every CI runner. That would leave the checked-in
R renderer and the completion gate untested in exactly the environment that is supposed to
prove they work, so these build the bundle by hand instead. The bundle is the renderer's only
input, so a hand-written one exercises the same path a real run does.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil

import pytest

from cmm.reporting import (
    render_transformation_report,
    validate_transformation_run,
)
from cmm.reporting.transformation import (
    TransformationReportError,
    renderer_script_path,
)

RANKING = """target_id,score,rank,bTS,mTS,wTS
PGI,3.5934,1,0.1369,0.1312,-0.1369
PFK,2.0970,2,0.1041,0.1007,-0.1040
TPI,1.9430,3,0.0975,0.0942,-0.0974
ENO,0.4210,4,0.0210,0.0203,-0.0210
"""
BASELINE = """target_id,moma_score,rank
ENO,-3.5072,1
PGI,-14.8486,2
PFK,-50.0562,3
TPI,-61.2200,4
"""
SWEEP = """epsilon,target_id,score,rank,rank_at_configured_epsilon
0.005,PGI,3.5934,1,1
0.005,PFK,2.0972,2,2
0.005,TPI,1.9431,3,3
0.02,PGI,3.5930,1,1
0.02,TPI,2.0100,2,3
0.02,PFK,1.9800,3,2
"""


def _r_is_ready() -> bool:
    """Whether this machine can run the checked-in renderer at all."""

    import subprocess

    if shutil.which("Rscript") is None:
        return False
    probe = (
        'q(status = as.integer(!all(vapply(c("jsonlite", "ggplot2", "ggrepel", '
        '"svglite", "ragg"), requireNamespace, logical(1), quietly = TRUE))))'
    )
    return subprocess.run(["Rscript", "--vanilla", "-e", probe]).returncode == 0


def _write(root: Path, relative: str, content: str) -> dict[str, object]:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    raw = path.read_bytes()
    return {
        "path": relative,
        "stage": relative.split("/")[0] if "/" in relative else "root",
        "media_type": "text/csv" if relative.endswith(".csv") else "application/json",
        "status": "complete",
        "sha256": hashlib.sha256(raw).hexdigest(),
        "size_bytes": len(raw),
    }


@pytest.fixture()
def bundle(tmp_path: Path) -> Path:
    """A minimal but complete transformation run, written directly to disk."""

    root = tmp_path / "run"
    root.mkdir()
    artifacts: dict[str, dict[str, object]] = {}

    def add(role: str, relative: str, content: str) -> None:
        entry = _write(root, relative, content)
        entry["role"] = role
        artifacts[role] = entry

    add("model", "model/model.xml", "<sbml/>\n")
    add(
        "preflight",
        "01_preflight/preflight.csv",
        "check,status,value,message\nsolver,pass,miqp,\n",
    )
    add(
        "source_reference_fluxes",
        "02_reference/source_reference_fluxes.csv",
        "reaction_id,flux\nPGI,4.86\nPFK,7.48\n",
    )
    add(
        "gene_differential_expression",
        "03_direction/gene_differential_expression.csv",
        "gene_id,log2_fold_change,t_statistic,p_value,significant,direction\nb4025,-3.0,-51.4,1e-6,True,-1\n",
    )
    add(
        "reaction_direction_map",
        "03_direction/reaction_direction_map.csv",
        "reaction_id,direction\nPGI,-1\nPFK,0\n",
    )
    add(
        "transformation_candidates",
        "04_candidates/candidates.csv",
        "target_id\nPGI\nPFK\nTPI\nENO\n",
    )
    add(
        "transformation_ranking",
        "05_transformation/transformation_ranking.csv",
        RANKING,
    )
    add("moma_baseline", "06_validation/moma_baseline.csv", BASELINE)
    add("epsilon_sensitivity", "06_validation/epsilon_sensitivity.csv", SWEEP)

    summary = {
        "method": "rmta",
        "n_candidates": 4,
        "top_target": "PGI",
        "top_score": 3.5934,
        "reference_method": "eflux2",
        "candidate_construction": {"n_open": 4, "n_blocked_removed": 0},
        "direction_construction": {"significance": "ttest", "n_reactions_changed": 20},
    }
    provenance = {
        "model_path": "model/model.xml",
        "model_id": "fixture",
        "solver": "gurobi",
        "method": "rmta",
        "perturbation": "reaction",
        "reference_method": "eflux2",
        "source_expression": "source.csv",
        "target_expression": "target.csv",
        "alpha": 0.66,
        "epsilon": 0.01,
        "reference_state_deviation": (
            "v_ref is a deterministic eflux2 solve; Yizhak et al. (2013) use iMAT, "
            "which CMM does not implement"
        ),
    }
    add("summary", "00_summary.json", json.dumps(summary, indent=2))
    add("provenance", "00_provenance.json", json.dumps(provenance, indent=2))
    add(
        "workflow_configuration",
        "00_config.json",
        json.dumps({"method": "rmta"}, indent=2),
    )

    (root / "00_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 2,
                "workflow": "transformation_target_discovery",
                "status": "complete",
                "artifacts": artifacts,
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return root


def test_a_bundle_validates_before_it_has_been_rendered(bundle: Path) -> None:
    """An analysis-only run is valid; it just says it has not been rendered."""

    report = validate_transformation_run(bundle)
    assert report.valid, report.issues
    assert report.phase == "pre-render"
    assert any("has not been rendered" in warning for warning in report.warnings)


@pytest.mark.skipif(
    not _r_is_ready(), reason="Rscript renderer packages are not installed"
)
def test_r_renderer_writes_vector_and_raster_for_every_panel(bundle: Path) -> None:
    report = render_transformation_report(bundle, highlight="PGI")

    manifest = json.loads(report.figure_manifest.read_text(encoding="utf-8"))
    rendered = {
        figure["id"] for figure in manifest["figures"] if figure["status"] == "rendered"
    }
    assert rendered == {
        "fig01_transformation_ranking",
        "fig02_ranking_vs_moma",
        "fig03_epsilon_sensitivity",
    }
    assert manifest["renderer"]["engine"] == "R/ggplot2"
    # The renderer records what produced the artwork, so a figure can be traced to a script.
    assert (
        manifest["renderer"]["script_sha256"]
        == hashlib.sha256(renderer_script_path().read_bytes()).hexdigest()
    )

    for figure in manifest["figures"]:
        for suffix in ("png", "pdf", "svg"):
            output = bundle / figure["outputs"][suffix]
            assert output.is_file() and output.stat().st_size > 0
        assert figure["dpi"] == 300

    linked = report.report_html.read_text(encoding="utf-8")
    standalone = report.report_standalone_html.read_text(encoding="utf-8")
    # The two variants exist precisely because one of them survives being moved.
    assert "src='figures/fig01_transformation_ranking.png'" in linked
    assert "src='figures/" not in standalone
    assert standalone.count("data:image/png;base64,") == 3

    assert validate_transformation_run(bundle).valid


@pytest.mark.skipif(
    not _r_is_ready(), reason="Rscript renderer packages are not installed"
)
def test_optional_panels_are_absent_with_a_stated_reason(
    bundle: Path, tmp_path: Path
) -> None:
    """A stage that was switched off must be legible as a choice, not read as a failure."""

    manifest_path = bundle / "00_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for role in ("moma_baseline", "epsilon_sensitivity"):
        (bundle / str(manifest["artifacts"][role]["path"])).write_text(
            "", encoding="utf-8"
        )
        manifest["artifacts"][role].update(
            {
                "status": "skipped",
                "reason": "the stage was disabled",
                "sha256": None,
                "size_bytes": 0,
            }
        )
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )

    report = render_transformation_report(bundle)
    figures = json.loads(report.figure_manifest.read_text(encoding="utf-8"))["figures"]
    unavailable = {f["id"]: f for f in figures if f["status"] != "rendered"}
    assert set(unavailable) == {"fig02_ranking_vs_moma", "fig03_epsilon_sensitivity"}
    for figure in unavailable.values():
        assert figure["reason"].strip()

    page = report.report_html.read_text(encoding="utf-8")
    assert "Not run" in page
    validation = validate_transformation_run(bundle)
    assert validation.valid, validation.issues


@pytest.mark.skipif(
    not _r_is_ready(), reason="Rscript renderer packages are not installed"
)
def test_a_required_panel_that_cannot_be_drawn_fails_the_render(bundle: Path) -> None:
    """The ranking panel is not optional, so an unusable ranking is an error, not a gap."""

    manifest_path = bundle / "00_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    path = bundle / str(manifest["artifacts"]["transformation_ranking"]["path"])
    path.write_text("target_id,score,rank\n", encoding="utf-8")
    entry = manifest["artifacts"]["transformation_ranking"]
    entry["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    entry["size_bytes"] = path.stat().st_size
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )

    with pytest.raises(TransformationReportError):
        render_transformation_report(bundle)


def test_render_refuses_a_directory_that_is_not_a_transformation_run(
    tmp_path: Path,
) -> None:
    with pytest.raises(TransformationReportError, match="no 00_manifest.json"):
        render_transformation_report(tmp_path)
