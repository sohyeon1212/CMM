"""Feature module boundaries and service entry points for the CMM workbench."""

from cmm.features.comparison import ComparisonResult, moma, reference_flux, room
from cmm.features.production import (
    fseof,
    fvseof,
    production_envelope,
    theoretical_yield,
)
from cmm.features.revert import revert_targets
from cmm.features.strain_design import (
    StrainDesign,
    StrainDesignResult,
    optknock,
    robustknock,
)
from cmm.features.transformation import direction_from_states, transformation_targets

INCLUDED_FEATURES = (
    "flux_visualization_slider",
    "fba",
    "fva",
    "dynamic_fba",
    "moma",
    "room",
    "batch_moma_room",
    "random_flux_sampling",
    "fseof",
    "fvseof",
    "omics_integration",
    "enzyme_constrained_modeling",
    "flux_response_analysis",
    "optknock",
    "robustknock",
    "revert_metabolism",
)

EXCLUDED_FEATURES = (
    "scenario_templates",
    "media_management",
    "scenario_file_formats",
    "elementary_flux_modes",
    "thermodynamics",
    "legacy_minimal_cut_sets",
)

__all__ = [
    "INCLUDED_FEATURES",
    "EXCLUDED_FEATURES",
    "ComparisonResult",
    "StrainDesign",
    "StrainDesignResult",
    "direction_from_states",
    "fseof",
    "fvseof",
    "moma",
    "production_envelope",
    "reference_flux",
    "optknock",
    "robustknock",
    "room",
    "revert_targets",
    "theoretical_yield",
    "transformation_targets",
]
