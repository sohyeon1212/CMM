from __future__ import annotations

import json
from dataclasses import replace

import pytest

import numpy as np
import pandas as pd
from cobra.io import write_sbml_model

from cmm.core.condition import Condition, ReactionBound
from cmm.workflows.transformation import (
    PUBLISHED_CHANGED_SET_RANGE,
    CandidateConfig,
    DirectionConfig,
    TransformationValidationConfig,
    TransformationWorkflowConfig,
    TransformationWorkflowError,
    _build_candidates,
    _check_expression_scale,
    _gene_directions,
    _read_expression,
    run_transformation_target_discovery,
)


def _config(**overrides) -> TransformationWorkflowConfig:
    values = dict(
        model_path="model.xml",
        source_expression_path="source.csv",
        target_expression_path="target.csv",
        imat_delta=0.5,
    )
    values.update(overrides)
    return TransformationWorkflowConfig(**values)


# --- defaults must be coherent with each other -------------------------------


def test_defaults_construct_and_follow_the_published_changed_set_size():
    config = _config()
    assert config.method == "mta"
    assert config.perturbation == "gene"
    assert config.top_n_changed == 200
    low, high = PUBLISHED_CHANGED_SET_RANGE
    assert low <= config.top_n_changed <= high
    assert config.follows_published_changed_set_size


def test_coupled_set_collapse_follows_the_perturbation_level():
    # Coupled sets are defined on reactions, so the default must not demand them of a
    # gene-level run — that combination used to make the default config unconstructible.
    assert _config().candidates.collapse_for("gene") is False
    assert _config(perturbation="reaction").candidates.collapse_for("reaction") is True
    # An explicit setting still wins in the direction that is coherent.
    forced_off = CandidateConfig(collapse_coupled_sets=False)
    assert forced_off.collapse_for("reaction") is False


def test_asking_for_coupled_sets_on_a_gene_run_is_rejected_with_the_way_out():
    with pytest.raises(ValueError, match="perturbation='reaction'"):
        _config(
            perturbation="gene",
            candidates=CandidateConfig(collapse_coupled_sets=True),
        )


# --- inputs ------------------------------------------------------------------


def test_source_and_target_must_differ():
    with pytest.raises(ValueError, match="different files"):
        _config(source_expression_path="same.csv", target_expression_path="same.csv")


@pytest.mark.parametrize(
    "overrides, message",
    [
        ({"alpha": 1.5}, "alpha"),
        ({"epsilon": -1.0}, "epsilon"),
        ({"parameter_k": 0.0}, "parameter_k"),
        ({"method": "rmta_continuous"}, "must be 'mta' or 'rmta'"),
        ({"perturbation": "enzyme"}, "must be 'gene' or 'reaction'"),
        ({"imat_delta": None}, "imat_delta is required"),
        ({"imat_delta": -0.1}, "imat_delta must be non-negative"),
        ({"imat_beta": 0.0}, "imat_beta"),
        ({"imat_beta": 1.5}, "imat_beta"),
        ({"imat_epsilon": 0.0}, "imat_epsilon"),
    ],
)
def test_invalid_values_are_rejected(overrides, message):
    with pytest.raises(ValueError, match=message):
        _config(**overrides)


@pytest.mark.parametrize(
    "key, value, message",
    [
        ("reference_method", "eflux2", "always built from iMAT"),
        ("reference_objective_fraction", 1.0, "no objective on growth"),
        ("epsilon_samples", 2000, "renamed to reference_samples"),
        ("epsilon_thinning", 500, "renamed to reference_thinning"),
    ],
)
def test_a_config_from_the_removed_eflux2_path_is_rejected_with_the_reason(
    key, value, message
):
    # An archived E-Flux2 config must not run as if it meant the iMAT reference.
    payload = {
        "model_path": "m.xml",
        "source_expression_path": "s.csv",
        "target_expression_path": "t.csv",
        "imat_delta": 0.5,
        key: value,
    }
    with pytest.raises(ValueError, match=message):
        TransformationWorkflowConfig.from_mapping(payload)


# --- direction section -------------------------------------------------------


def test_p_value_ranking_requires_the_t_test():
    with pytest.raises(ValueError, match="ranking='p_value' needs"):
        _config(
            direction=DirectionConfig(significance="fold_change", ranking="p_value")
        )


def test_fold_change_run_may_rank_on_fold_change():
    config = _config(
        direction=DirectionConfig(significance="fold_change", ranking="fold_change")
    )
    assert config.direction.significance == "fold_change"


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"p_value_cutoff": 0.0}, "p_value_cutoff"),
        ({"p_value_cutoff": 1.5}, "p_value_cutoff"),
        ({"up_threshold": -1.0}, "fold-change"),
        ({"top_n_changed": 0}, "top_n_changed"),
    ],
)
def test_direction_section_validates(kwargs, message):
    with pytest.raises(ValueError, match=message):
        _config(direction=DirectionConfig(**kwargs))


def test_no_cut_is_allowed_but_leaves_the_published_range():
    config = _config(direction=DirectionConfig(top_n_changed=None))
    assert config.top_n_changed is None
    assert not config.follows_published_changed_set_size


# --- epsilon suggestion ------------------------------------------------------


def test_suggest_epsilon_reads_percentiles_off_the_reference_state():
    fluxes = {f"r{i}": float(i) for i in range(1, 101)}
    suggestion = TransformationWorkflowConfig.suggest_epsilon(fluxes)
    assert suggestion["p10"] < suggestion["median"] < suggestion["p75"]
    assert suggestion["median"] == pytest.approx(50.5, abs=1.0)


def test_suggest_epsilon_ignores_zero_flux_and_uses_magnitude():
    # Sign is irrelevant to "how far must this move"; a zero carries no scale information.
    assert TransformationWorkflowConfig.suggest_epsilon({"a": 0.0, "b": 0.0}) == {}
    signed = TransformationWorkflowConfig.suggest_epsilon({"a": -4.0, "b": 4.0})
    assert signed["median"] == pytest.approx(4.0)


# --- provenance --------------------------------------------------------------


def test_provenance_records_how_v_ref_was_built_and_what_still_differs():
    provenance = _config(imat_delta=0.7).to_provenance()
    assert provenance["reference_method"] == "imat"
    assert provenance["imat_delta"] == 0.7
    assert provenance["imat_beta"] == 1.0
    assert provenance["imat_epsilon"] == 1.0
    assert provenance["imat_gpr_rule"] == "AND=min, OR=max"
    assert provenance["reference_samples"] == 2000
    assert provenance["reference_thinning"] == 500
    deviations = provenance["reference_state_deviations"]
    assert len(deviations) == 3
    joined = " ".join(deviations)
    assert "AND=min, OR=max" in joined
    assert "warm-up" in joined and "one seeded chain" in joined
    # Using one of iMAT's optima is the papers' own procedure, so it is not listed.
    assert "optim" not in joined


def test_provenance_records_the_resolved_candidate_construction():
    # The candidate count is the denominator of any "top N%" claim, so how it was built rides
    # with the numbers rather than being inferable from the method name.
    gene = _config().to_provenance()
    reaction = _config(perturbation="reaction").to_provenance()
    assert gene["candidate_collapse_coupled_sets"] is False
    assert reaction["candidate_collapse_coupled_sets"] is True


def test_provenance_carries_every_parameter_that_changes_a_ranking():
    provenance = _config(
        alpha=0.5, epsilon=0.01, method="rmta", perturbation="reaction"
    ).to_provenance()
    for key in ("alpha", "epsilon", "method", "perturbation", "top_n_changed"):
        assert key in provenance
    assert provenance["method"] == "rmta"
    assert provenance["epsilon"] == 0.01


# --- serialization -----------------------------------------------------------


def test_from_json_resolves_relative_paths_against_the_config_file(tmp_path):
    (tmp_path / "sub").mkdir()
    payload = {
        "model_path": "sub/model.xml",
        "source_expression_path": "sub/source.csv",
        "target_expression_path": "sub/target.csv",
        "imat_delta": 0.5,
        "output_dir": "runs/out",
        "method": "rmta",
        "epsilon": 0.01,
        "direction": {"significance": "ttest", "top_n_changed": 150},
        "candidates": {"essential_growth_fraction": 0.1},
        "validation": {"epsilon_sweep": [0.001, 0.01]},
    }
    path = tmp_path / "config.json"
    path.write_text(json.dumps(payload), encoding="utf-8")

    config = TransformationWorkflowConfig.from_json(path)
    assert config.model_path == (tmp_path / "sub" / "model.xml").resolve()
    assert config.output_dir == (tmp_path / "runs" / "out").resolve()
    assert config.method == "rmta"
    assert config.direction.top_n_changed == 150
    assert config.candidates.essential_growth_fraction == 0.1
    assert config.validation.epsilon_sweep == (0.001, 0.01)


def test_from_json_rejects_a_non_object(tmp_path):
    path = tmp_path / "config.json"
    path.write_text("[1, 2]", encoding="utf-8")
    with pytest.raises(ValueError, match="must contain an object"):
        TransformationWorkflowConfig.from_json(path)


def test_from_mapping_builds_a_condition():
    config = TransformationWorkflowConfig.from_mapping(
        {
            "model_path": "m.xml",
            "source_expression_path": "s.csv",
            "target_expression_path": "t.csv",
            "imat_delta": 0.5,
            "condition": {
                "name": "anaerobic",
                "bounds": [{"reaction_id": "EX_o2_e", "lower_bound": 0.0}],
            },
        }
    )
    assert isinstance(config.condition, Condition)
    assert config.condition.name == "anaerobic"
    assert config.condition.bounds == (
        ReactionBound(reaction_id="EX_o2_e", lower_bound=0.0),
    )


def test_validation_section_rejects_a_negative_sweep_value():
    with pytest.raises(ValueError, match="epsilon_sweep"):
        _config(validation=TransformationValidationConfig(epsilon_sweep=(-0.1,)))


# --- execution ---------------------------------------------------------------


def _write_expression(path, genes, values, columns=("r1", "r2", "r3")):
    pd.DataFrame(values, index=genes, columns=list(columns)).to_csv(path)


@pytest.fixture
def transformation_inputs(tmp_path, branched_model):
    """Real six-reaction solves fit the restricted CI license, including MOMA and MIQP."""

    model_path = tmp_path / "model.xml"
    write_sbml_model(branched_model, str(model_path))
    genes = [gene.id for gene in branched_model.genes]
    rng = np.random.default_rng(0)
    levels = {"g1": 50.0, "g2": 100.0, "g3": 1.0, "g5": 1.0, "gb": 50.0}
    source = np.array([levels[gene] for gene in genes])[:, None] * rng.uniform(
        0.95, 1.05, (len(genes), 3)
    )
    target = source.copy()
    target[genes.index("g2")] *= 0.01
    target[genes.index("g3")] *= 100.0
    target[genes.index("g5")] *= 100.0
    _write_expression(tmp_path / "source.csv", genes, source)
    _write_expression(tmp_path / "target.csv", genes, target)
    return model_path, tmp_path / "source.csv", tmp_path / "target.csv"


def test_a_non_sbml_model_is_refused_before_anything_runs(tmp_path, branched_model):
    from cobra.io import save_json_model

    json_model = tmp_path / "model.json"
    save_json_model(branched_model, str(json_model))
    config = _config(model_path=json_model)
    with pytest.raises(ValueError, match=r"reads SBML models.*load_json_model"):
        run_transformation_target_discovery(config)
    with pytest.raises(FileNotFoundError, match="not a file"):
        run_transformation_target_discovery(_config(model_path=tmp_path / "absent.xml"))


def test_read_expression_rejects_a_table_with_no_numeric_columns(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text("gene,label\nb0001,up\n", encoding="utf-8")
    with pytest.raises(TransformationWorkflowError, match="no numeric"):
        _read_expression(path)


def test_read_expression_rejects_duplicate_gene_ids(tmp_path):
    path = tmp_path / "dup.csv"
    path.write_text("gene,a\nb0001,1\nb0001,2\n", encoding="utf-8")
    with pytest.raises(TransformationWorkflowError, match="repeats a gene id"):
        _read_expression(path)


def test_t_test_needs_replicates_and_says_what_to_do_instead(tmp_path):
    genes = ["g1", "g2"]
    _write_expression(tmp_path / "s.csv", genes, np.ones((2, 1)), columns=("only",))
    _write_expression(tmp_path / "t.csv", genes, np.zeros((2, 1)), columns=("only",))
    source = _read_expression(tmp_path / "s.csv")
    target = _read_expression(tmp_path / "t.csv")
    # The workflow re-raises cmm.omics' message as its own error type, so the caller sees
    # both what went wrong and the supported way round it.
    with pytest.raises(TransformationWorkflowError, match=r"gene_directions\(\)"):
        _gene_directions(source, target, DirectionConfig(significance="ttest"))


def test_gene_directions_are_signed_toward_the_target():
    source = pd.DataFrame(np.full((2, 3), 8.0), index=["up", "down"])
    target = source.copy()
    target.loc["up"] += 3.0
    target.loc["down"] -= 3.0
    target += np.array([[0.01, -0.01, 0.0], [0.01, -0.01, 0.0]])
    frame = _gene_directions(source, target, DirectionConfig())
    assert frame.loc["up", "direction"] == 1
    assert frame.loc["down", "direction"] == -1


def test_disjoint_gene_identifiers_are_a_stop_not_a_silent_empty_result():
    source = pd.DataFrame(np.ones((2, 3)), index=["a", "b"])
    target = pd.DataFrame(np.ones((2, 3)), index=["x", "y"])
    with pytest.raises(TransformationWorkflowError, match="share no gene ids"):
        _gene_directions(source, target, DirectionConfig())


def test_linear_fold_change_is_a_ratio_not_a_difference():
    source = pd.DataFrame(
        {"measurement": [100.0, 0.0, 3.0]}, index=["small", "on", "off"]
    )
    target = pd.DataFrame({"measurement": [102.0, 3.0, 0.0]}, index=source.index)
    config = DirectionConfig(significance="fold_change", ranking="fold_change")
    evidence = _gene_directions(source, target, config)
    assert evidence.loc["small", "log2_fold_change"] == pytest.approx(
        np.log2(103 / 101)
    )
    assert evidence["direction"].to_dict() == {"small": 0, "on": 1, "off": -1}
    assert evidence.loc["on", "log2_fold_change"] == 2
    assert evidence.loc["off", "log2_fold_change"] == -2
    assert source.loc["small", "measurement"] == 100.0


def test_ttest_uses_log2_of_each_linear_replicate():
    from cmm.omics.differential import gene_directions_from_replicates

    source = pd.DataFrame([[1.0, 5.0, 20.0]], index=["g1"])
    target = pd.DataFrame([[2.0, 40.0, 120.0]], index=["g1"])
    expected = gene_directions_from_replicates(np.log2(source + 1), np.log2(target + 1))
    pd.testing.assert_frame_equal(
        _gene_directions(source, target, DirectionConfig()), expected
    )


@pytest.mark.parametrize("value", [float("inf"), float("-inf")])
def test_expression_inputs_require_finite_measurements(tmp_path, value):
    path = tmp_path / "expression.csv"
    pd.DataFrame({"measurement": [value]}, index=["g1"]).to_csv(path)
    for scale in ("linear", "log2"):
        with pytest.raises(TransformationWorkflowError, match="finite"):
            _read_expression(path, scale)


def test_negative_expression_is_rejected_as_linear_but_accepted_as_log2(tmp_path):
    path = tmp_path / "expression.csv"
    pd.DataFrame({"measurement": [-1.0]}, index=["g1"]).to_csv(path)
    with pytest.raises(TransformationWorkflowError, match="expression_scale='log2'"):
        _read_expression(path)
    assert _read_expression(path, "log2").loc["g1", "measurement"] == -1.0


def test_candidate_filters_preserve_distinct_reaction_and_gene_effects(
    parallel_pathway_model,
):
    model = parallel_pathway_model
    reaction_candidates, _ = _build_candidates(model, _config(perturbation="reaction"))
    assert reaction_candidates == ("R2", "R3")
    gene_candidates, record = _build_candidates(model, _config())
    assert gene_candidates == ("g2", "g3")
    assert record["n_genes_essential_removed"] == 1
    assert model.slim_optimize() == pytest.approx(10)
    assert model.reactions.BIOMASS.lower_bound == 1

    unfiltered, _ = _build_candidates(
        model, _config(candidates=CandidateConfig(exclude_essential=False))
    )
    assert unfiltered == ("g1", "g2", "g3")


def test_gene_essentiality_checks_joint_deletion_of_nonessential_reactions(
    parallel_pathway_model,
):
    model = parallel_pathway_model
    model.reactions.R2.gene_reaction_rule = "g2 and joint"
    model.reactions.R3.gene_reaction_rule = "g3 and joint"
    candidates, record = _build_candidates(model, _config())
    assert candidates == ("g2", "g3")
    assert record["n_genes_essential_removed"] == 2
    assert model.slim_optimize() == pytest.approx(10)


@pytest.fixture
def small_transformation_config(tmp_path, parallel_pathway_model):
    model_path = tmp_path / "model.xml"
    write_sbml_model(parallel_pathway_model, str(model_path))
    # Identical basenames must not cause the archived source and target to overwrite each other.
    paths = []
    for state, value in (("source", 100.0), ("target", 1.0)):
        directory = tmp_path / state
        directory.mkdir()
        path = directory / "expression.tsv"
        pd.DataFrame({"measurement": [value] * 3}, index=["g1", "g2", "g3"]).to_csv(
            path, sep="\t"
        )
        paths.append(path)
    return _config(
        model_path=model_path,
        source_expression_path=paths[0],
        target_expression_path=paths[1],
        output_dir=tmp_path / "run",
        epsilon=0.01,
        reference_thinning=10,  # sampling speed only; the default is 500
        candidates=CandidateConfig(explicit=("g2",)),
        direction=DirectionConfig(
            significance="fold_change", ranking="fold_change", top_n_changed=None
        ),
        validation=TransformationValidationConfig(enabled=False),
    )


@pytest.mark.requires_miqp
def test_workflow_does_not_promote_a_two_percent_linear_change(
    small_transformation_config,
):
    config = small_transformation_config
    pd.DataFrame({"measurement": [102.0] * 3}, index=["g1", "g2", "g3"]).to_csv(
        config.target_expression_path, sep="\t"
    )
    with pytest.raises(
        TransformationWorkflowError, match="no reaction was labelled as changed"
    ):
        run_transformation_target_discovery(config)


@pytest.mark.requires_miqp
def test_bundle_replays_after_relocation_without_original_inputs(
    small_transformation_config, tmp_path
):
    import shutil
    import subprocess
    import sys
    from pathlib import Path

    from cmm.reporting import validate_transformation_run

    config = small_transformation_config
    result = run_transformation_target_discovery(config)
    root = result.run_directory
    manifest = json.loads((root / "00_manifest.json").read_text())
    for role, source in (
        ("model", config.model_path),
        ("source_expression", config.source_expression_path),
        ("target_expression", config.target_expression_path),
    ):
        assert (
            root / manifest["artifacts"][role]["path"]
        ).read_bytes() == source.read_bytes()
        source.unlink()

    relocated = tmp_path / "relocated"
    shutil.move(root, relocated)
    assert validate_transformation_run(relocated).valid
    reproduced = tmp_path / "reproduced"
    process = subprocess.run(
        [
            sys.executable,
            str(relocated / "scripts/reproduce.py"),
            "--output-dir",
            str(reproduced),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert process.returncode == 0, process.stderr
    assert validate_transformation_run(reproduced).valid
    ranking_path = "05_transformation/transformation_ranking.csv"
    assert (reproduced / ranking_path).read_bytes() == (
        relocated / ranking_path
    ).read_bytes()

    # A rerun whose inputs are its own archived copies must read them before overwrite cleanup.
    archived = TransformationWorkflowConfig.from_json(relocated / "00_config.json")
    before = {
        field: Path(getattr(archived, field)).read_bytes()
        for field in ("model_path", "source_expression_path", "target_expression_path")
    }
    run_transformation_target_discovery(replace(archived, overwrite=True))
    for field, data in before.items():
        assert Path(getattr(archived, field)).read_bytes() == data
    assert validate_transformation_run(relocated).valid


@pytest.mark.requires_miqp
def test_workflow_baseline_preserves_failed_moma_status(small_transformation_config):
    config = replace(
        small_transformation_config,
        candidates=CandidateConfig(explicit=("g1", "g2")),
        validation=TransformationValidationConfig(run_moma_baseline=True),
    )
    result = run_transformation_target_discovery(config)
    baseline = pd.read_csv(
        result.run_directory / "06_validation/moma_baseline.csv"
    ).set_index("target_id")
    assert baseline.loc["g1", "status"] == "infeasible"
    assert baseline.loc["g1", "moma_score"] == float("-inf")
    assert baseline.loc["g2", "status"] == "optimal"


@pytest.mark.requires_miqp
def test_run_writes_a_schema_v2_bundle_with_one_role_per_artifact(
    tmp_path, transformation_inputs
):
    model_path, source, target = transformation_inputs
    output = tmp_path / "run"
    config = TransformationWorkflowConfig(
        imat_delta=0.5,
        reference_thinning=10,  # sampling speed only; the default is 500
        model_path=model_path,
        source_expression_path=source,
        target_expression_path=target,
        output_dir=output,
        perturbation="reaction",
        epsilon=0.01,
        direction=DirectionConfig(top_n_changed=20),
        validation=TransformationValidationConfig(epsilon_sweep=(0.001,)),
    )
    result = run_transformation_target_discovery(config)

    assert result.run_directory == output.resolve()
    manifest = json.loads((output / "00_manifest.json").read_text())
    assert manifest["schema_version"] == 2
    assert manifest["workflow"] == "transformation_target_discovery"
    for role in (
        "model",
        "preflight",
        "source_reference_fluxes",
        "gene_differential_expression",
        "reaction_direction_map",
        "transformation_candidates",
        "transformation_ranking",
        "moma_baseline",
        "epsilon_sensitivity",
        "provenance",
        "summary",
        "workflow_configuration",
    ):
        assert role in manifest["artifacts"], role
        assert (output / manifest["artifacts"][role]["path"]).is_file()

    ranking = pd.read_csv(output / "05_transformation/transformation_ranking.csv")
    assert len(ranking) == len(result.candidates)
    assert list(ranking["rank"]) == sorted(ranking["rank"])
    assert ranking["score"].is_monotonic_decreasing


@pytest.mark.requires_miqp
def test_run_records_how_the_candidate_set_was_built(tmp_path, transformation_inputs):
    # The candidate count is the denominator of any "top N%" reading, so the construction
    # must be recoverable from the run rather than inferred from the method name.
    model_path, source, target = transformation_inputs
    config = TransformationWorkflowConfig(
        imat_delta=0.5,
        reference_thinning=10,  # sampling speed only; the default is 500
        model_path=model_path,
        source_expression_path=source,
        target_expression_path=target,
        perturbation="reaction",
        epsilon=0.01,
        direction=DirectionConfig(top_n_changed=20),
        validation=TransformationValidationConfig(enabled=False),
    )
    result = run_transformation_target_discovery(config)
    built = result.candidate_filtering
    assert built["source"] == "constructed"
    assert built["n_reactions_allowed"] >= len(result.candidates)
    assert built["coupling"]["coupling"] == "full"
    assert result.summary()["n_candidates"] == len(result.candidates)


@pytest.mark.requires_miqp
def test_explicit_candidates_skip_construction(tmp_path, transformation_inputs):
    model_path, source, target = transformation_inputs
    config = TransformationWorkflowConfig(
        imat_delta=0.5,
        reference_thinning=10,  # sampling speed only; the default is 500
        model_path=model_path,
        source_expression_path=source,
        target_expression_path=target,
        perturbation="reaction",
        epsilon=0.01,
        direction=DirectionConfig(top_n_changed=20),
        candidates=CandidateConfig(explicit=("R2", "R3", "R5")),
        validation=TransformationValidationConfig(enabled=False),
    )
    result = run_transformation_target_discovery(config)
    assert set(result.candidates) == {"R2", "R3", "R5"}
    assert result.candidate_filtering["source"] == "explicit"


@pytest.mark.requires_miqp
def test_refusing_to_overwrite_a_non_empty_directory(tmp_path, transformation_inputs):
    model_path, source, target = transformation_inputs
    output = tmp_path / "run"
    output.mkdir()
    (output / "someone_elses_file.txt").write_text("keep me", encoding="utf-8")
    config = TransformationWorkflowConfig(
        imat_delta=0.5,
        reference_thinning=10,  # sampling speed only; the default is 500
        model_path=model_path,
        source_expression_path=source,
        target_expression_path=target,
        output_dir=output,
        perturbation="reaction",
        candidates=CandidateConfig(explicit=("R2",)),
        validation=TransformationValidationConfig(enabled=False),
    )
    with pytest.raises(FileExistsError, match="not empty"):
        run_transformation_target_discovery(config)
    assert (output / "someone_elses_file.txt").read_text() == "keep me"


# --- report rendering --------------------------------------------------------


@pytest.mark.requires_miqp
def test_report_renders_figures_and_states_what_it_must(
    tmp_path, transformation_inputs
):
    from cmm.reporting import (
        render_transformation_report,
        validate_transformation_run,
    )

    model_path, source, target = transformation_inputs
    output = tmp_path / "run"
    config = TransformationWorkflowConfig(
        imat_delta=0.5,
        reference_thinning=10,  # sampling speed only; the default is 500
        model_path=model_path,
        source_expression_path=source,
        target_expression_path=target,
        output_dir=output,
        method="rmta",
        perturbation="reaction",
        epsilon=0.01,
        direction=DirectionConfig(top_n_changed=20),
        candidates=CandidateConfig(explicit=("R2", "R3", "R5")),
        validation=TransformationValidationConfig(epsilon_sweep=(0.001,)),
    )
    result = run_transformation_target_discovery(config)
    report = render_transformation_report(result.run_directory, highlight="R2")

    assert report.report_html.is_file()
    names = {path.name for path in report.figures}
    assert "fig02_transformation_ranking.png" in names
    assert "fig03_ranking_vs_moma.png" in names
    assert "fig04_epsilon_sensitivity.png" in names
    for path in report.figures:
        assert path.stat().st_size > 0
        # Editable vector output beside the raster: a figure that exists only as a PNG has to
        # be redrawn before it can go anywhere else.
        for suffix in (".svg", ".pdf"):
            assert path.with_suffix(suffix).stat().st_size > 0

    # The linked page references figures relatively, so it renders blank once it is moved --
    # and says nothing about it. The standalone copy is the one that survives being sent.
    standalone = report.report_standalone_html.read_text(encoding="utf-8")
    assert "data:image/png;base64," in standalone
    assert "src='figures/" not in standalone

    # A rendered page is not a finished run; the gate is what says so.
    validation = validate_transformation_run(result.run_directory)
    assert validation.valid, validation.issues
    assert validation.phase == "post-render"

    page = report.report_html.read_text(encoding="utf-8")
    # A reader must not have to open the provenance file to learn any of these.
    assert "iMAT" in page
    assert "How the reference state was built" in page
    assert "OptGP samples, thinning 10" in page  # the run's own sampling settings
    assert "warm-up" in page  # what still differs from the papers is stated
    assert "not a finding" in page  # the source/target direction is an input
    assert "denominator" in page  # the candidate count qualifies any percentile
    assert "chosen, not derived" in page  # epsilon
    assert "in silico" in page
    assert "Yizhak" in page and "Valc" in page  # both methods cited
    assert "R2" in page
    # The decomposition replaces the default baseline: same question, no extra solve.
    assert "Score decomposition" in page
    assert "exchanges permit uptake" in page
    assert "model genes with an expression value" in page
    assert "not comparable across runs" in page


@pytest.mark.requires_miqp
def test_rmta_run_decomposes_its_score_without_a_baseline(
    tmp_path, transformation_inputs
):
    model_path, source, target = transformation_inputs
    config = TransformationWorkflowConfig(
        imat_delta=0.5,
        reference_thinning=10,  # sampling speed only; the default is 500
        model_path=model_path,
        source_expression_path=source,
        target_expression_path=target,
        output_dir=tmp_path / "run",
        method="rmta",
        perturbation="reaction",
        epsilon=0.01,
        direction=DirectionConfig(top_n_changed=20),
        candidates=CandidateConfig(explicit=("R2", "R3", "R5")),
    )
    assert config.validation.run_moma_baseline is False
    result = run_transformation_target_discovery(config)
    assert result.moma_baseline == ()

    by_id = {row["target_id"]: row for row in result.ranking}
    decomposition = {row["target_id"]: row for row in result.score_decomposition}
    assert set(decomposition) == set(by_id)
    for target_id, row in decomposition.items():
        ranked = by_id[target_id]
        assert row["rTS"] == ranked["score"]
        assert row["rank_rTS"] == ranked["rank"]
        assert row["amplified"] == (
            row["bTS"] > 0 and row["mTS"] > 0 and row["wTS"] < 0
        )
        assert row["rank_delta"] == row["rank_mTS"] - row["rank_rTS"]
        # Equation 9: outside the amplification branch the score is the MOMA component.
        if not row["amplified"]:
            assert row["rTS"] == row["mTS"]
    assert sorted(row["rank_mTS"] for row in decomposition.values()) == [1, 2, 3]

    structure = result.summary()["ranking_structure"]
    assert structure["n_amplified"] == sum(
        1 for row in decomposition.values() if row["amplified"]
    )
    assert {"n_distinct_scores", "largest_tie_block", "spearman_rTS_vs_mTS"} <= set(
        structure
    )
    coverage = result.summary()["data_coverage"]
    assert coverage["genes_measured_in_model"] <= coverage["model_genes"]
    assert 0 <= coverage["n_exchanges_open_for_uptake"] <= coverage["n_exchanges"]

    manifest = json.loads((result.run_directory / "00_manifest.json").read_text())
    assert manifest["artifacts"]["score_decomposition"]["status"] == "complete"
    assert manifest["artifacts"]["moma_baseline"]["status"] == "skipped"


@pytest.mark.requires_miqp
def test_mta_run_skips_the_decomposition_with_a_reason(tmp_path, transformation_inputs):
    from cmm.core.media import Medium

    model_path, source, target = transformation_inputs
    config = TransformationWorkflowConfig(
        imat_delta=0.5,
        reference_thinning=10,  # sampling speed only; the default is 500
        model_path=model_path,
        source_expression_path=source,
        target_expression_path=target,
        output_dir=tmp_path / "run",
        method="mta",
        perturbation="reaction",
        epsilon=0.01,
        direction=DirectionConfig(top_n_changed=20),
        candidates=CandidateConfig(explicit=("R2", "R3")),
        medium=Medium(name="custom", uptake={"SUP_A": 10.0}),
    )
    result = run_transformation_target_discovery(config)
    assert result.score_decomposition == ()
    assert "n_amplified" not in result.ranking_structure
    entry = json.loads((result.run_directory / "00_manifest.json").read_text())[
        "artifacts"
    ]["score_decomposition"]
    assert entry["status"] == "skipped" and "only for rMTA" in entry["reason"]
    # A Medium object used to be recorded as None, which the report printed as "model as
    # loaded" although the run was constrained.
    provenance = json.loads((result.run_directory / "00_provenance.json").read_text())
    assert provenance["medium"] == "custom"
    assert provenance["medium_applied"]["applied"] == {"SUP_A": 10.0}


@pytest.mark.requires_miqp
def test_mta_run_does_not_publish_three_copies_of_one_score(
    tmp_path, transformation_inputs
):
    # CMM returns the same value in all four score slots for method="mta". Emitting bTS/mTS/wTS
    # would present one number as three independent measurements that happen to agree.
    model_path, source, target = transformation_inputs
    config = TransformationWorkflowConfig(
        imat_delta=0.5,
        reference_thinning=10,  # sampling speed only; the default is 500
        model_path=model_path,
        source_expression_path=source,
        target_expression_path=target,
        method="mta",
        perturbation="reaction",
        epsilon=0.01,
        direction=DirectionConfig(top_n_changed=20),
        candidates=CandidateConfig(explicit=("R2", "R3")),
        validation=TransformationValidationConfig(enabled=False),
    )
    rows = run_transformation_target_discovery(config).ranking
    assert not {"bTS", "mTS", "wTS"} & set(rows[0])
    # The display name is unrelated to the component question and rides along wherever the
    # model names the thing being knocked out, so the guard is on the components, not on the
    # exact column set.
    assert set(rows[0]) <= {"target_id", "target_name", "score", "rank"}


@pytest.mark.requires_miqp
def test_rmta_run_publishes_the_three_components(tmp_path, transformation_inputs):
    # Equation 9 branches on their signs, so a reader cannot reconstruct which branch fired
    # from the combined score alone.
    model_path, source, target = transformation_inputs
    config = TransformationWorkflowConfig(
        imat_delta=0.5,
        reference_thinning=10,  # sampling speed only; the default is 500
        model_path=model_path,
        source_expression_path=source,
        target_expression_path=target,
        method="rmta",
        perturbation="reaction",
        epsilon=0.01,
        direction=DirectionConfig(top_n_changed=20),
        candidates=CandidateConfig(explicit=("R2", "R3")),
        validation=TransformationValidationConfig(enabled=False),
    )
    rows = run_transformation_target_discovery(config).ranking
    assert {"bTS", "mTS", "wTS"} <= set(rows[0])


def test_report_refuses_a_run_from_a_different_workflow(tmp_path):
    from cmm.reporting import TransformationReportError, render_transformation_report

    (tmp_path / "00_manifest.json").write_text(
        json.dumps({"workflow": "production_target_discovery", "artifacts": {}}),
        encoding="utf-8",
    )
    with pytest.raises(TransformationReportError, match="not a transformation run"):
        render_transformation_report(tmp_path)


def test_report_refuses_a_directory_with_no_manifest(tmp_path):
    from cmm.reporting import TransformationReportError, render_transformation_report

    with pytest.raises(TransformationReportError, match="no 00_manifest.json"):
        render_transformation_report(tmp_path)


@pytest.mark.requires_miqp
def test_completion_gate_catches_the_failures_a_browser_shows_as_success(
    tmp_path, transformation_inputs
):
    """Each of these renders as a page that opens, which is why the gate exists."""

    import re
    import shutil

    from cmm.reporting import render_transformation_report, validate_transformation_run

    model_path, source, target = transformation_inputs
    config = TransformationWorkflowConfig(
        imat_delta=0.5,
        reference_thinning=10,  # sampling speed only; the default is 500
        model_path=model_path,
        source_expression_path=source,
        target_expression_path=target,
        output_dir=tmp_path / "run",
        method="mta",
        perturbation="reaction",
        epsilon=0.01,
        direction=DirectionConfig(top_n_changed=20),
        candidates=CandidateConfig(explicit=("R2", "R3", "R5")),
    )
    result = run_transformation_target_discovery(config)
    render_transformation_report(result.run_directory)
    assert validate_transformation_run(result.run_directory).valid

    def mutated(name: str, mutate) -> tuple[bool, tuple[str, ...]]:
        copy = tmp_path / name
        shutil.rmtree(copy, ignore_errors=True)
        shutil.copytree(result.run_directory, copy)
        mutate(copy)
        report = validate_transformation_run(copy)
        return report.valid, report.issues

    # A CSV edited after the run: every number the report quotes is now unattributable.
    valid, issues = mutated(
        "edited",
        lambda root: (root / "05_transformation/transformation_ranking.csv").write_text(
            (root / "05_transformation/transformation_ranking.csv")
            .read_text()
            .replace("R2", "XXX"),
            encoding="utf-8",
        ),
    )
    assert not valid and any("sha256" in issue for issue in issues)

    # A figure the manifest still claims was drawn. The page shows an empty box.
    valid, issues = mutated(
        "no_png",
        lambda root: (root / "figures/fig02_transformation_ranking.png").unlink(),
    )
    assert not valid and any("missing non-empty png" in issue for issue in issues)

    # The standalone copy pointing at a file that will not travel with it.
    valid, issues = mutated(
        "delinked",
        lambda root: (root / "report_standalone.html").write_text(
            re.sub(
                r"src='data:image/png;base64,[^']*'",
                "src='figures/fig02_transformation_ranking.png'",
                (root / "report_standalone.html").read_text(encoding="utf-8"),
                count=1,
            ),
            encoding="utf-8",
        ),
    )
    assert not valid and any("does not carry" in issue for issue in issues)


def test_report_refuses_a_production_run_directory(tmp_path):
    from cmm.reporting import validate_transformation_run

    run = tmp_path / "run"
    run.mkdir()
    (run / "00_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 2,
                "workflow": "production_target_discovery",
                "artifacts": {},
            }
        ),
        encoding="utf-8",
    )
    report = validate_transformation_run(run)
    assert not report.valid
    assert any("not a transformation_target_discovery run" in i for i in report.issues)


# --- per-reaction epsilon, expression scale, FDR ------------------------------


def test_epsilon_defaults_to_the_published_per_reaction_rule():
    config = _config()
    assert config.epsilon == "per_reaction"
    assert config.reference_samples == 2000
    provenance = config.to_provenance()
    assert provenance["epsilon"] == "per_reaction"
    # v_ref needs the samples either way, so a scalar run records them too.
    assert _config(epsilon=0.01).to_provenance()["reference_samples"] == 2000


@pytest.mark.parametrize(
    "overrides, message",
    [
        ({"epsilon": "auto"}, "a number or 'per_reaction'"),
        ({"reference_samples": 1}, "reference_samples"),
        ({"reference_thinning": 0}, "reference_thinning"),
        ({"expression_scale": "ln"}, "expression_scale"),
        (
            {"direction": DirectionConfig(p_adjust="bonferroni")},  # type: ignore[arg-type]
            "p_adjust must be",
        ),
        (
            {
                "direction": DirectionConfig(
                    significance="fold_change", ranking="fold_change", p_adjust="fdr_bh"
                )
            },
            "significance='ttest' only",
        ),
    ],
)
def test_epsilon_scale_and_adjustment_settings_are_validated(overrides, message):
    with pytest.raises(ValueError, match=message):
        _config(**overrides)


def test_mapping_round_trip_keeps_the_new_fields():
    config = TransformationWorkflowConfig.from_mapping(
        {
            "model_path": "m.xml",
            "source_expression_path": "s.csv",
            "target_expression_path": "t.csv",
            "imat_delta": 0.5,
            "expression_scale": "log2",
            "epsilon": "per_reaction",
            "reference_thinning": 1000,
            "direction": {"p_adjust": "fdr_bh"},
        }
    )
    assert config.expression_scale == "log2"
    assert config.reference_thinning == 1000
    assert config.direction.p_adjust == "fdr_bh"
    provenance = config.to_provenance()
    assert provenance["p_adjust"] == "fdr_bh"
    assert provenance["imat_discretisation_transform"] == "value"
    assert provenance["direction_transform"] == "value"


@pytest.mark.parametrize(
    "scale, values, ok, phrase",
    [
        ("linear", [0.0, 50.0, 1200.0], True, "consistent with linear"),
        ("linear", [2.0, 8.5, 15.9], False, "set expression_scale='log2'"),
        ("log2", [-1.0, 8.5, 15.9], True, "consistent with log2"),
        ("log2", [0.0, 50.0, 1200.0], False, "set expression_scale='linear'"),
    ],
)
def test_preflight_flags_a_file_whose_values_contradict_its_declared_scale(
    scale, values, ok, phrase
):
    records = []

    def check(name, passed, value, message):
        records.append((name, passed, value, message))

    table = pd.DataFrame({"r1": values, "r2": values}, index=["a", "b", "c"])
    _check_expression_scale(check, "source", table, scale)
    ((name, passed, value, message),) = records
    assert name == "source_expression_scale"
    assert passed is ok
    assert value == {"declared": scale, "min": min(values), "max": max(values)}
    assert phrase in message


def test_a_log2_file_gives_the_same_directions_as_its_linear_original():
    rng = np.random.default_rng(3)
    genes = [f"g{i}" for i in range(12)]
    source = pd.DataFrame(rng.uniform(10, 100, (12, 3)), index=genes)
    target = source * rng.uniform(0.2, 5.0, (12, 1))
    from_linear = _gene_directions(source, target, DirectionConfig())
    from_log2 = _gene_directions(
        np.log2(source + 1.0), np.log2(target + 1.0), DirectionConfig(), "log2"
    )
    pd.testing.assert_frame_equal(from_linear, from_log2)


@pytest.mark.requires_miqp
def test_per_reaction_run_writes_thresholds_bounds_and_the_range_check(
    tmp_path, transformation_inputs
):
    from cmm.reporting import validate_transformation_run

    model_path, source, target = transformation_inputs
    output = tmp_path / "run"
    config = TransformationWorkflowConfig(
        imat_delta=0.5,
        reference_thinning=10,  # sampling speed only; the default is 500
        model_path=model_path,
        source_expression_path=source,
        target_expression_path=target,
        output_dir=output,
        perturbation="reaction",
        reference_samples=200,
        direction=DirectionConfig(top_n_changed=20),
        validation=TransformationValidationConfig(enabled=False),
    )
    result = run_transformation_target_discovery(config)

    changed = {rid for rid, value in result.direction.items() if value != 0}
    thresholds = pd.read_csv(output / "05_transformation/reaction_epsilons.csv")
    assert set(thresholds["reaction_id"]) == changed
    assert (thresholds["epsilon"] >= 1e-3).all()
    assert len(result.reaction_epsilons) == len(changed)

    # v_ref is the mean of the same 200 samples the thresholds came from.
    assert result.provenance["reference_state"]["parameters"]["n_samples"] == 200
    pattern = pd.read_csv(output / "02_reference/imat_reference_pattern.csv")
    assert set(pattern["reaction_id"]) == set(result.reference_fluxes)
    assert dict(zip(pattern["reaction_id"], pattern["v_ref"])) == pytest.approx(
        dict(result.reference_fluxes)
    )
    states = pd.read_csv(output / "02_reference/gene_expression_states.csv")
    assert set(states["state"]) <= {-1, 0, 1}

    ranking_provenance = result.provenance["ranking_provenance"]
    assert ranking_provenance["epsilon"] == "per_reaction"
    assert result.summary()["direction_construction"]["epsilon_summary"][
        "n_reactions"
    ] == len(changed)

    # Six reactions cannot reach the published 100-200, and the run says so.
    size_check = next(r for r in result.preflight if r.check == "changed_set_size")
    assert size_check.status == "warning"
    assert size_check.value == len(changed)
    assert (
        result.summary()["direction_construction"]["changed_set_in_published_range"]
        is False
    )

    manifest = json.loads((output / "00_manifest.json").read_text())
    assert manifest["artifacts"]["reaction_epsilons"]["status"] == "complete"
    assert validate_transformation_run(output).phase  # schema accepted


@pytest.mark.requires_miqp
def test_scalar_run_marks_the_threshold_table_skipped(tmp_path, transformation_inputs):
    model_path, source, target = transformation_inputs
    output = tmp_path / "run"
    config = TransformationWorkflowConfig(
        imat_delta=0.5,
        reference_thinning=10,  # sampling speed only; the default is 500
        model_path=model_path,
        source_expression_path=source,
        target_expression_path=target,
        output_dir=output,
        perturbation="reaction",
        epsilon=0.01,
        reference_samples=200,
        direction=DirectionConfig(top_n_changed=20),
        validation=TransformationValidationConfig(enabled=False),
    )
    result = run_transformation_target_discovery(config)
    manifest = json.loads((output / "00_manifest.json").read_text())
    assert manifest["artifacts"]["reaction_epsilons"]["status"] == "skipped"
    # The reference is sampled whatever epsilon is, so its pattern is always written.
    assert manifest["artifacts"]["imat_reference_pattern"]["status"] == "complete"
    assert result.reaction_epsilons == ()


def test_minimum_fold_change_is_validated_and_recorded():
    config = _config(direction=DirectionConfig(min_abs_log2fc=0.263))
    assert config.to_provenance()["min_abs_log2fc"] == 0.263
    with pytest.raises(ValueError, match="non-negative"):
        _config(direction=DirectionConfig(min_abs_log2fc=-0.1))
    with pytest.raises(ValueError, match="up_threshold/down_threshold"):
        _config(
            direction=DirectionConfig(
                significance="fold_change", ranking="fold_change", min_abs_log2fc=0.5
            )
        )


def test_reader_reports_which_columns_became_samples(tmp_path):
    path = tmp_path / "expression.csv"
    pd.DataFrame(
        {
            "symbol": ["a", "b", "c"],
            "entrez": [101, 202, 303],
            "r1": [1.5, 2.5, 2.5],
            "r2": [1.6, 2.4, 2.6],
        },
        index=["g1", "g2", "g3"],
    ).to_csv(path)
    frame = _read_expression(path)
    assert list(frame.columns) == ["entrez", "r1", "r2"]
    columns = frame.attrs["columns"]
    assert columns.dropped == ("symbol",)
    assert columns.identifier_like == ("entrez",)


@pytest.mark.requires_miqp
def test_preflight_lists_sample_columns_and_flags_identifier_columns(
    tmp_path, transformation_inputs
):
    model_path, source, target = transformation_inputs
    with_id = pd.read_csv(source, index_col=0)
    with_id.insert(0, "entrez", range(1000, 1000 + len(with_id)))
    flagged = tmp_path / "source_with_id.csv"
    with_id.to_csv(flagged)
    result = run_transformation_target_discovery(
        TransformationWorkflowConfig(
            imat_delta=0.5,
            reference_thinning=10,  # sampling speed only; the default is 500
            model_path=model_path,
            source_expression_path=flagged,
            target_expression_path=target,
            perturbation="reaction",
            epsilon=0.01,
            direction=DirectionConfig(top_n_changed=20, min_abs_log2fc=0.5),
            validation=TransformationValidationConfig(enabled=False),
        )
    )
    records = {record.check: record for record in result.preflight}
    assert records["target_sample_columns"].status == "pass"
    assert records["target_sample_columns"].value["samples"] == ["r1", "r2", "r3"]
    source_record = records["source_sample_columns"]
    assert source_record.status == "warning"
    assert source_record.value["identifier_like"] == ["entrez"]
    assert "identifiers" in source_record.message
    assert result.summary()["direction_construction"]["min_abs_log2fc"] == 0.5


@pytest.mark.requires_miqp
def test_a_log2_file_gives_the_reference_and_moma_baseline_of_its_linear_original(
    tmp_path, transformation_inputs
):
    # Every place expression is used must put it on one scale. A linear file is discretised as
    # log2(x + 1), so a file holding exactly those values must give the same iMAT pattern, the
    # same seeded samples, and so the same v_ref and MOMA-baseline target.
    model_path, source, target = transformation_inputs
    log2_paths = {}
    for name, path in (("source", source), ("target", target)):
        frame = pd.read_csv(path, index_col=0)
        log2_paths[name] = tmp_path / f"{name}_log2.csv"
        np.log2(frame + 1.0).to_csv(log2_paths[name])

    def run(src, tgt, scale):
        return run_transformation_target_discovery(
            TransformationWorkflowConfig(
                imat_delta=0.5,
                reference_thinning=10,  # sampling speed only; the default is 500
                model_path=model_path,
                source_expression_path=src,
                target_expression_path=tgt,
                expression_scale=scale,
                perturbation="reaction",
                epsilon=0.01,
                direction=DirectionConfig(top_n_changed=20),
                validation=TransformationValidationConfig(run_moma_baseline=True),
            )
        )

    linear = run(source, target, "linear")
    logged = run(log2_paths["source"], log2_paths["target"], "log2")
    for result in (linear, logged):
        scale_checks = [r for r in result.preflight if r.check.endswith("_scale")]
        assert len(scale_checks) == 2
        assert all(r.status == "pass" for r in scale_checks)
    for rid, flux in linear.reference_fluxes.items():
        assert logged.reference_fluxes[rid] == pytest.approx(flux, abs=1e-6)
    assert logged.moma_baseline, "the MOMA baseline must run on log2 input"
    by_id = {row["target_id"]: row["moma_score"] for row in linear.moma_baseline}
    for row in logged.moma_baseline:
        assert row["moma_score"] == pytest.approx(by_id[row["target_id"]], abs=1e-6)
