"""Validated, deterministic publication-report artifacts for CMM scenario runs."""

from cmm.reporting.publication import (
    FigureManifest,
    PublicationBundle,
    ReportBuildResult,
    build_publication_bundle,
    build_publication_report,
    render_publication_figures,
    render_production_report,
)
from cmm.reporting.transformation import (
    TransformationReport,
    TransformationReportError,
    render_transformation_report,
)
from cmm.reporting.schema import (
    ARTIFACT_CONTRACTS,
    RUN_SCHEMA_VERSION,
    RunArtifact,
    RunValidationError,
    ValidationReport,
    ValidatedRun,
    validate_production_run,
    validate_run,
)

__all__ = [
    "ARTIFACT_CONTRACTS",
    "TransformationReport",
    "TransformationReportError",
    "render_transformation_report",
    "RUN_SCHEMA_VERSION",
    "FigureManifest",
    "PublicationBundle",
    "ReportBuildResult",
    "RunArtifact",
    "RunValidationError",
    "ValidationReport",
    "ValidatedRun",
    "build_publication_bundle",
    "build_publication_report",
    "render_publication_figures",
    "render_production_report",
    "validate_production_run",
    "validate_run",
]
