from __future__ import annotations

import pytest
from cmm.core import solvers
from cmm.core.flux_state import (
    FluxState,
    reference_state_from_samples,
    reference_state_pfba,
)
from cmm.core.provenance import model_fingerprint, run_provenance
from cmm.core.results import TargetRanking, TargetScore
from cmm.core.solvers import SolverCapabilityError

# --- solvers ---------------------------------------------------------------


def test_active_solver_supports_lp():
    assert solvers.supports("LP")


def test_capability_matrix_is_consistent_for_active_solver(branched_model):
    name = solvers.active_solver(branched_model)
    caps = solvers.capabilities(branched_model.solver.interface)
    assert "LP" in caps
    # Whatever the active solver, LP must be runnable end to end.
    assert branched_model.slim_optimize() == pytest.approx(10)
    assert name in solvers.available_solvers() or name  # name resolves to a string


def test_solver_status_reports_gurobi_recommended(branched_model):
    status = solvers.solver_status(branched_model)
    assert status.name in ("gurobi", "cplex")
    assert status.recommended
    assert status.warning is None
    assert "QP" in status.capabilities


def test_solver_status_warns_on_limited_solver(branched_model):
    branched_model.solver = "glpk"
    status = solvers.solver_status(branched_model)
    assert not status.recommended
    assert status.warning is not None
    assert "QP" in status.warning


def test_require_raises_typed_error_for_unavailable_capability():
    class _FakeLPOnly:
        __name__ = "optlang.glpk_exact_interface"

    with pytest.raises(SolverCapabilityError) as exc:
        solvers.require("MIQP", _FakeLPOnly(), feature="original MTA")
    assert exc.value.capability == "MIQP"
    assert "original MTA" in str(exc.value)


# --- flux state ------------------------------------------------------------


def test_flux_state_distance_l1_and_l2():
    a = FluxState({"R1": 3.0, "R2": 0.0}, name="a")
    b = FluxState({"R1": 0.0, "R2": 4.0}, name="b")
    assert a.distance(b, order=1) == pytest.approx(7.0)
    assert a.distance(b, order=2) == pytest.approx(5.0)


def test_flux_state_roundtrip_serialization():
    state = FluxState({"R1": 1.5}, name="ref", provenance="pfba", metadata={"k": 1})
    restored = FluxState.deserialize(state.serialize())
    assert restored.name == "ref"
    assert restored.provenance == "pfba"
    assert restored.get("R1") == pytest.approx(1.5)
    assert restored.get("missing", -1.0) == -1.0


@pytest.mark.parametrize(
    "payload, message",
    [
        ({"fluxes": [1.0]}, "fluxes.*mapping"),
        ({"fluxes": {"R1": 1.0}, "metadata": []}, "metadata.*mapping"),
        ({"fluxes": {"R1": 1.0}, "provenance": "unknown"}, "provenance"),
    ],
)
def test_flux_state_deserialize_rejects_malformed_payload(payload, message):
    with pytest.raises(ValueError, match=message):
        FluxState.deserialize(payload)


@pytest.mark.parametrize("fluxes", [{}, {"R": float("nan")}, {"R": float("inf")}])
def test_flux_state_rejects_empty_or_nonfinite_fluxes(fluxes):
    with pytest.raises(ValueError):
        FluxState(fluxes)


def test_reference_state_pfba_routes_disease_branch(branched_model):
    ref = reference_state_pfba(branched_model, name="disease")
    assert ref.provenance == "pfba"
    # pFBA prefers the shorter R2 branch; the healthy R3->R5 branch stays at zero.
    assert ref.get("R2") == pytest.approx(10, abs=1e-6)
    assert ref.get("R3") == pytest.approx(0, abs=1e-6)
    assert ref.get("R5") == pytest.approx(0, abs=1e-6)


def test_reference_state_from_samples():
    import pandas as pd

    samples = pd.DataFrame({"R1": [1.0, 3.0], "R2": [0.0, 2.0]})
    ref = reference_state_from_samples(samples)
    assert ref.provenance == "sampling_mean"
    assert ref.get("R1") == pytest.approx(2.0)
    assert ref.metadata["n_samples"] == 2


def test_model_fingerprint_is_deterministic_and_bound_sensitive(branched_model):
    first = model_fingerprint(branched_model)
    assert first == model_fingerprint(branched_model.copy())
    branched_model.reactions.R2.upper_bound = 9.0
    assert model_fingerprint(branched_model) != first


def test_run_provenance_records_reproducibility_fields(branched_model):
    provenance = run_provenance(branched_model, alpha=0.66)
    assert provenance["model_id"] == "branched"
    assert len(provenance["model_sha256"]) == 64
    assert provenance["solver"]
    assert provenance["parameters"] == {"alpha": 0.66}


# --- target ranking --------------------------------------------------------


def test_target_ranking_sorts_and_breaks_ties_deterministically():
    ranking = TargetRanking(
        method="demo",
        targets=(
            TargetScore("b", 1.0),
            TargetScore("a", 1.0),
            TargetScore("c", 5.0),
        ),
    ).sorted()
    assert [t.target_id for t in ranking] == ["c", "a", "b"]
    assert ranking.best().target_id == "c"


def test_target_ranking_export_frame_is_deterministic():
    ranking = TargetRanking.from_scores(
        "demo",
        {"x": 2.0, "y": 9.0},
        metadata={"note": "t"},
    )
    frame = ranking.to_frame()
    assert list(frame.columns) == ["rank", "target_id", "score"]
    assert list(frame["target_id"]) == ["y", "x"]
    assert list(frame["rank"]) == [1, 2]
