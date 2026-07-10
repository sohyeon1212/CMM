"""Non-optional validation on the 2,583-reaction E. coli iJO1366 reconstruction."""

from __future__ import annotations

import cobra.flux_analysis as cfa
import pytest
from cobra.io import load_model

from cmm.core import fba, fva, pfba
from cmm.features.production import fseof, theoretical_yield


@pytest.fixture(scope="module")
def ijo1366():
    model = load_model("iJO1366")
    # The open LP backend has no size-limited model license and makes this CI-reproducible.
    model.solver = "glpk"
    return model


@pytest.mark.genome_scale
def test_ijo1366_is_a_real_genome_scale_model(ijo1366):
    assert len(ijo1366.reactions) == 2583
    assert len(ijo1366.genes) == 1367


@pytest.mark.genome_scale
def test_genome_scale_fba_matches_direct_cobra(ijo1366):
    result = fba(ijo1366)
    direct = ijo1366.optimize()
    assert result.status == "optimal"
    assert result.objective_value == pytest.approx(direct.objective_value, abs=1e-8)
    assert result.metadata["model_sha256"]


@pytest.mark.genome_scale
def test_genome_scale_pfba_matches_direct_cobra(ijo1366):
    result = pfba(ijo1366)
    direct = cfa.pfba(ijo1366)
    assert result.status == "optimal"
    assert sum(abs(value) for value in result.fluxes.values()) == pytest.approx(
        direct.fluxes.abs().sum(), rel=1e-7
    )


@pytest.mark.genome_scale
def test_genome_scale_fva_matches_direct_cobra(ijo1366):
    reaction_ids = ["BIOMASS_Ec_iJO1366_core_53p95M", "EX_glc__D_e", "EX_succ_e"]
    result = fva(ijo1366, reactions=reaction_ids, fraction_of_optimum=0.9)
    direct = cfa.flux_variability_analysis(
        ijo1366, reaction_list=reaction_ids, fraction_of_optimum=0.9, processes=1
    )
    for rid in reaction_ids:
        assert result[rid].minimum == pytest.approx(
            direct.loc[rid, "minimum"], abs=1e-7
        )
        assert result[rid].maximum == pytest.approx(
            direct.loc[rid, "maximum"], abs=1e-7
        )


@pytest.mark.genome_scale
def test_genome_scale_production_scan_returns_actionable_targets(ijo1366):
    product = "EX_succ_e"
    biomass = "BIOMASS_Ec_iJO1366_core_53p95M"
    yield_result = theoretical_yield(ijo1366, product, substrate="EX_glc__D_e")
    scan = fseof(
        ijo1366,
        product,
        biomass,
        n_steps=3,
        reactions=["PPC", "FRD2", product, biomass],
    )
    assert yield_result.status == "optimal"
    assert scan.metadata["n_failed_levels"] == 0
    assert product not in scan.amplification_targets()
    assert biomass not in scan.amplification_targets()
