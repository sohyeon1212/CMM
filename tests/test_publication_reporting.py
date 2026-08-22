from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

import cmm.reporting.publication as publication

from cmm.reporting import (
    RunValidationError,
    build_publication_report,
    render_production_report,
    render_publication_figures,
    validate_production_run,
    validate_run,
)
from cmm.reporting.publication import (
    FigureRenderError,
    _amplification_support_statement,
    _decode_renderer_stream,
    _renderer_environment,
    renderer_script_path,
)


def _write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _record(path: Path, root: Path, role: str) -> dict[str, object]:
    media_type = (
        "text/csv"
        if path.suffix == ".csv"
        else "application/json"
        if path.suffix == ".json"
        else "text/x-python"
    )
    return {
        "path": path.relative_to(root).as_posix(),
        "stage": path.parent.name,
        "role": role,
        "media_type": media_type,
        "status": "complete",
        "sha256": _digest(path),
        "size_bytes": path.stat().st_size,
    }


def _entry(
    root: Path,
    relative: str,
    role: str,
    supplementary: list[dict[str, object]],
) -> dict[str, object]:
    path = root / relative
    entry = {
        "path": relative,
        "status": "complete",
        "sha256": _digest(path),
        "size_bytes": path.stat().st_size,
    }
    if path.suffix == ".csv":
        metadata_path = path.with_suffix(".metadata.json")
        _write(
            metadata_path,
            json.dumps(
                {
                    "model_sha256": "a" * 64,
                    "method": role,
                    "parameters": {"fixture": True},
                },
                sort_keys=True,
            )
            + "\n",
        )
        entry["metadata_path"] = metadata_path.relative_to(root).as_posix()
        supplementary.append(_record(metadata_path, root, f"{role}_metadata"))
    return entry


def _replace_artifact(
    root: Path,
    manifest: dict[str, object],
    role: str,
    content: str,
) -> None:
    artifacts = manifest["artifacts"]
    assert isinstance(artifacts, dict)
    entry = artifacts[role]
    assert isinstance(entry, dict)
    path = root / str(entry["path"])
    _write(path, content)
    entry["sha256"] = _digest(path)
    entry["size_bytes"] = path.stat().st_size


def _update_csv_artifact_rows(
    root: Path,
    manifest: dict[str, object],
    role: str,
    *,
    target_column: str,
    target: str,
    updates: dict[str, str],
) -> None:
    artifacts = manifest["artifacts"]
    assert isinstance(artifacts, dict)
    entry = artifacts[role]
    assert isinstance(entry, dict)
    path = root / str(entry["path"])
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert rows
    fields = list(rows[0])
    matched = False
    for row in rows:
        if row[target_column] == target:
            row.update(updates)
            matched = True
    assert matched
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    _replace_artifact(root, manifest, role, output.getvalue())


def _make_run(
    root: Path,
    *,
    optional: str = "complete",
    empty_design: bool = False,
    many_targets: bool = False,
    exhaustive_validation: bool = False,
) -> Path:
    provenance = {
        "model_id": "fixture_model",
        "model_sha256": "a" * 64,
        "source_model_sha256": "b" * 64,
        "conditioned_model_sha256": "a" * 64,
        "solver": "fixture-solver",
        "solver_version": "1.0",
        "python": "3.12",
        "cmm": "0.test",
        "cobra": "0.test",
        "numpy": "0.test",
        "pandas": "0.test",
        "scipy": "0.test",
        "parameters": {
            "medium": "defined_medium",
            "condition": "anaerobic_fixture",
            "product": "EX_product",
            "reference_method": "pfba",
        },
        "medium_application": {
            "name": "defined_medium",
            "mode": "explicit",
            "applied": {"EX_substrate": 10.0, "EX_o2": 0.0},
            "dropped": [],
        },
    }
    summary = {
        "schema_version": 2,
        "status": "complete",
        "model_id": "fixture_model",
        "model_sha256": "a" * 64,
        "product": "EX_product",
        "substrate": "EX_substrate",
        "biomass": "BIOMASS",
        "wild_type_growth": 0.4,
        "wild_type_product": 0.1,
        "theoretical_product_flux": 12.0,
        "theoretical_molar_yield": 1.2,
        "warnings": [
            "summary oxygen-bound warning",
            "The model medium was used as loaded; review exchange bounds.",
        ],
    }
    config = {
        "model_path": "model/fixture.xml",
        "product": "EX_product",
        "substrate": "EX_substrate",
        "reference_method": "pfba",
        "viability_fraction": 0.1,
        "envelope_points": 3,
        "fseof_steps": 3,
        "fvseof_steps": 3,
        "validation": {"sampling": {"n": 2, "method": "achr", "seed": 7}},
    }
    if exhaustive_validation:
        amplification_count = 20 if many_targets else 2
        response_count = amplification_count + 9
        coverage = {
            "single_knockout_candidates_expected": 9,
            "amplification_candidates_expected": amplification_count,
            "flux_response_expected": response_count,
            "flux_response_attempted": response_count,
            "flux_response_completed": response_count,
            "flux_response_failed": 0,
            "sampling_expected": 10,
            "sampling_attempted": 10,
            "sampling_completed": 10,
            "sampling_failed": 0,
            "sampling_skipped": 0,
        }
        summary["validation_candidate_policy"] = {
            "single_gene_knockout": "all_display_ranked_candidates",
            "amplification": "all_report_selected_candidates",
        }
        summary["validation_coverage"] = coverage
        provenance["validation_coverage"] = coverage
    _write(root / "00_provenance.json", json.dumps(provenance, sort_keys=True) + "\n")
    _write(root / "00_summary.json", json.dumps(summary, sort_keys=True) + "\n")
    invocation_config = {
        **config,
        "model_path": str((root / "source-model.xml").resolve()),
        "output_dir": str(root.resolve()),
    }
    _write(
        root / "00_config.json", json.dumps(invocation_config, sort_keys=True) + "\n"
    )
    _write(root / "model" / "fixture.xml", "<sbml><model id='fixture'/></sbml>\n")
    _write(
        root / "01_preflight" / "preflight.csv",
        "check,status,value,message\n"
        "model_structure,pass,4 reactions,model structure is usable\n"
        "condition,warning,anaerobic,fixture condition warning\n"
        "medium,warning,model_as_loaded,model bounds were retained under the model-as-loaded medium\n",
    )
    fluxes = "reaction_id,flux\nBIOMASS,0.4\nEX_product,0.1\nR1,1.0\nR2,2.0\n"
    _write(root / "03_reference" / "wild_type_fluxes.csv", fluxes)
    _write(root / "03_reference" / "reference_pfba.csv", fluxes)
    _write(
        root / "02_yield" / "theoretical_yield.csv",
        "status,molar_yield,product_flux,substrate_uptake,carbon_imbalance\n"
        "optimal,1.2,12.0,10.0,false\n",
    )
    _write(
        root / "02_yield" / "production_envelope.csv",
        "product_flux,growth_min,growth_max\n0,0,0.4\n6,0,0.25\n12,0,0\n",
    )
    screen_header = (
        "method,method_rank,display_rank,selected,target_id,kind,blocked_reactions,"
        "blocked_reaction_signature,status,growth_rate,objective,growth_fraction,"
        "target_production,product_flux,product_delta,product_fold_change,improves_product,"
        "objective_value,distance,distance_kind,n_changed_reactions,n_blocked_reactions,"
        "target_name,blocked_reaction_names,blocked_reaction_equations,blocked_reaction_gprs\n"
    )
    _write(
        root / "04_single_knockout" / "moma.csv",
        screen_header
        + "moma_l2,2,2,true,g1,gene,R1,R1,optimal,0.31,0.31,0.775,3.2,3.2,3.1,32,"
        "true,2,1.4,euclidean,1,1,Gene 1,Reaction 1,A --> B,g1\n"
        + "moma_l2,1,1,true,g2,gene,R2,R2,optimal,0.27,0.27,0.675,4.1,4.1,4,41,"
        "true,3,1.7,euclidean,1,1,Gene 2,Reaction 2,B --> C,g2 and g4\n"
        + "moma_l2,,,false,g4,gene,R2,R2,optimal,0.28,0.28,0.7,4,4,3.9,40,true,"
        "2.9,1.6,euclidean,1,1,Gene 4,Reaction 2,B --> C,g2 and g4\n"
        + "moma_l2,,3,false,g5,gene,R5,R5,optimal,0.26,0.26,0.65,3,3,2.9,30,true,"
        "2.5,1.5,euclidean,1,1,Gene 5,Reaction 5,C --> D,g5\n"
        + "moma_l2,,4,false,g6,gene,R6,R6,optimal,0.24,0.24,0.6,2.5,2.5,2.4,25,true,"
        "2.1,1.3,euclidean,1,1,Gene 6,Reaction 6,D --> E,g6\n"
        + "moma_l2,,5,false,g7,gene,R7,R7,optimal,0.22,0.22,0.55,2,2,1.9,20,true,"
        "1.8,1.1,euclidean,1,1,Gene 7,Reaction 7,E --> F,g7\n"
        + "moma_l2,,,false,g3,gene,,,infeasible,,,,,,,,,,,,0,Gene 3,,,,\n",
    )
    _write(
        root / "04_single_knockout" / "room.csv",
        screen_header
        + "room,1,1,true,g1,gene,R1,R1,optimal,0.29,0.29,0.725,4.2,4.2,4.1,42,true,"
        "4,,,1,1,Gene 1,Reaction 1,A --> B,g1\n"
        + "room,2,,true,g2,gene,R2,R2,optimal,0.25,0.25,0.625,3.7,3.7,3.6,37,true,"
        "5,,,1,1,Gene 2,Reaction 2,B --> C,g2 and g4\n"
        + "room,,2,false,g4,gene,R2,R2,optimal,0.28,0.28,0.7,4,4,3.9,40,true,"
        "4.5,,,1,1,Gene 4,Reaction 2,B --> C,g2 and g4\n"
        + "room,,3,false,g5,gene,R5,R5,optimal,0.27,0.27,0.675,3.6,3.6,3.5,36,true,"
        "4.2,,,1,1,Gene 5,Reaction 5,C --> D,g5\n"
        + "room,,4,false,g6,gene,R6,R6,optimal,0.26,0.26,0.65,3.4,3.4,3.3,34,true,"
        "4.1,,,1,1,Gene 6,Reaction 6,D --> E,g6\n"
        + "room,,5,false,g7,gene,R7,R7,optimal,0.2,0.2,0.5,1.5,1.5,1.4,15,true,"
        "3.8,,,1,1,Gene 7,Reaction 7,E --> F,g7\n"
        + "room,,,false,g3,gene,,,infeasible,,,,,,,,,,,,0,Gene 3,,,,\n",
    )
    _write(
        root / "04_single_knockout" / "consensus.csv",
        "target_id,recommended,selected_methods\n"
        "g1,true,moma_l2;room\ng2,true,moma_l2;room\ng3,false,\n",
    )
    _write(
        root / "04_single_knockout" / "gene_knockout_mapping.csv",
        "gene_id,gene_name,inert,blocked_reaction,reaction_name,reaction_equation,gpr\n"
        "g1,Gene 1,false,R1,Reaction 1,A --> B,g1\n"
        "g2,Gene 2,false,R2,Reaction 2,B --> C,g2 and g4\n"
        "g4,Gene 4,false,R2,Reaction 2,B --> C,g2 and g4\n"
        "g5,Gene 5,false,R5,Reaction 5,C --> D,g5\n"
        "g6,Gene 6,false,R6,Reaction 6,D --> E,g6\n"
        "g7,Gene 7,false,R7,Reaction 7,E --> F,g7\n"
        "g_cross,Gene Cross,false,R2,Reaction 2,B --> C,g2 and g_cross\n"
        "g_cross,Gene Cross,false,RX,Reaction X,X --> Y,g_cross\n"
        "g3,Gene 3,true,,,,\n",
    )
    if exhaustive_validation:

        def screen_row(
            method: str, target: str, signature: str, display_rank: int | None
        ) -> str:
            values = [
                method,
                "",
                "" if display_rank is None else str(display_rank),
                "false",
                target,
                "gene",
                signature,
                signature,
                "optimal",
                "0.2",
                "0.2",
                "0.5",
                "1.0",
                "1.0",
                "0.9",
                "10",
                "true",
                "1.0",
                "0.5" if method == "moma_l2" else "",
                "euclidean" if method == "moma_l2" else "",
                "1",
                "1",
                f"Gene {target[1:]}",
                f"Reaction {signature}",
                f"{signature}_a --> {signature}_b",
                target,
            ]
            return ",".join(values) + "\n"

        moma_path = root / "04_single_knockout" / "moma.csv"
        room_path = root / "04_single_knockout" / "room.csv"
        moma_text = moma_path.read_text(encoding="utf-8") + "".join(
            screen_row("moma_l2", f"g{target}", f"R{target}", None)
            for target in range(8, 12)
        )
        room_text = room_path.read_text(encoding="utf-8")
        for rank, target in ((2, "g4"), (3, "g5"), (4, "g6"), (5, "g7")):
            room_text = room_text.replace(
                f"room,,{rank},false,{target}", f"room,,,false,{target}"
            )
        room_text += "".join(
            screen_row("room", f"g{target}", f"R{target}", target - 6)
            for target in range(8, 12)
        )
        _write(moma_path, moma_text)
        _write(room_path, room_text)
        mapping_path = root / "04_single_knockout" / "gene_knockout_mapping.csv"
        mapping_text = mapping_path.read_text(encoding="utf-8")
        mapping_text += (
            "g1_alias,Gene 1 alias,false,R1,Reaction 1,A --> B,g1 or g1_alias\n"
        )
        mapping_text += "".join(
            f"g{target},Gene {target},false,R{target},Reaction R{target},"
            f"R{target}_a --> R{target}_b,g{target}\n"
            for target in range(8, 12)
        )
        _write(mapping_path, mapping_text)
    design_header = "knockouts,growth,max_product,guaranteed_product,growth_coupled\n"
    if empty_design:
        _write(root / "05_strain_design" / "optknock.csv", design_header)
        _write(root / "05_strain_design" / "robustknock.csv", design_header)
    else:
        _write(
            root / "05_strain_design" / "optknock.csv",
            design_header + '"R1;R2",0.20,8.0,3.0,true\nR3,0.30,6.0,0.0,false\n',
        )
        _write(
            root / "05_strain_design" / "robustknock.csv",
            design_header + '"R1;R2",0.20,7.0,6.5,true\nR4,0.22,5.5,5.0,true\n',
        )

    if many_targets:
        fseof_rows = "".join(
            f"F{target},{step},{1000 + step if target == 1 else target + step},{1 / target},"
            f"amplification,true,fseof,{target},true,complete,"
            f"{'true' if target == 1 else 'false'},"
            f"{0.01 if target == 1 else 0.9},"
            f"{'loop-dominated diagnostic' if target == 1 else ''}\n"
            for target in range(1, 11)
            for step in (0, 6, 12)
        )
        fvseof_rows = "".join(
            f"V{target},{step},{1000 + step if target == 1 else target + step},"
            f"{900 + step if target == 1 else target + step - 0.5},10,amplification,1,true,"
            f"true,fvseof,{target},true,complete,{'true' if target == 1 else 'false'},"
            f"{0.01 if target == 1 else 0.9},"
            f"{'loop-dominated diagnostic' if target == 1 else ''}\n"
            for target in range(1, 11)
            for step in (0, 6, 12)
        )
        fseof_rank_rows = "".join(
            f"F{target},{target},{1 / target},amplify,true,fseof,{target},true,complete,"
            f"{'true' if target == 1 else 'false'},"
            f"{0.01 if target == 1 else 0.9},"
            f"{'loop-dominated diagnostic' if target == 1 else ''}\n"
            for target in range(1, 11)
        )
        fvseof_rank_rows = "".join(
            f"V{target},{target},amplify,1,true,{1 / target},-0.1,10,true,fvseof,{target},"
            f"true,complete,{'true' if target == 1 else 'false'},"
            f"{0.01 if target == 1 else 0.9},"
            f"{'loop-dominated diagnostic' if target == 1 else ''}\n"
            for target in range(1, 11)
        )
    else:
        fseof_rows = (
            "R1,0,1,0.4,amplification,true,fseof,1,true,complete,false,0.9,\n"
            "R1,6,3,0.4,amplification,true,fseof,1,true,complete,false,0.9,\n"
            "R1,12,5,0.4,amplification,true,fseof,1,true,complete,false,0.9,\n"
            "R2,0,2,0.3,amplification,true,fseof,2,true,complete,false,0.9,\n"
            "R2,6,3,0.3,amplification,true,fseof,2,true,complete,false,0.9,\n"
            "R2,12,4,0.3,amplification,true,fseof,2,true,complete,false,0.9,\n"
        )
        fvseof_rows = (
            "R1,0,1,0.5,4.5,amplification,1,true,true,fvseof,1,true,complete,false,0.9,\n"
            "R1,6,3,2,4.5,amplification,1,true,true,fvseof,1,true,complete,false,0.9,\n"
            "R1,12,5,4.5,4.5,amplification,1,true,true,fvseof,1,true,complete,false,0.9,\n"
            "R2,0,2,1,3.5,amplification,2,true,true,fvseof,2,true,complete,false,0.9,\n"
            "R2,6,3,2.5,3.5,amplification,2,true,true,fvseof,2,true,complete,false,0.9,\n"
            "R2,12,4,3.5,3.5,amplification,2,true,true,fvseof,2,true,complete,false,0.9,\n"
        )
        fseof_rank_rows = (
            "R1,1,0.4,amplify,true,fseof,1,true,complete,false,0.9,\n"
            "R2,2,0.3,amplify,true,fseof,2,true,complete,false,0.9,\n"
        )
        fvseof_rank_rows = (
            "R1,1,amplify,1,true,0.4,-0.1,5,true,fvseof,1,true,complete,false,0.9,\n"
            "R2,2,amplify,2,true,0.3,-0.1,4,true,fvseof,2,true,complete,false,0.9,\n"
        )
    _write(
        root / "06_amplification" / "fseof.csv",
        "reaction_id,amplification_rank,slope,classification,actionable,proposal_method,"
        "method_rank,report_selected,loop_diagnostic_status,loop_artifact_flag,"
        "loopless_to_standard_capacity_ratio,loop_diagnostic_reason\n"
        + fseof_rank_rows,
    )
    _write(
        root / "06_amplification" / "fvseof.csv",
        "reaction_id,amplification_rank,classification,park_type,robust,slope,capacity_slope,"
        "mean_capacity,actionable,proposal_method,method_rank,report_selected,"
        "loop_diagnostic_status,loop_artifact_flag,loopless_to_standard_capacity_ratio,"
        "loop_diagnostic_reason\n" + fvseof_rank_rows,
    )
    _write(
        root / "06_amplification" / "fseof_tidy.csv",
        "target,enforced_product_flux,reaction_flux,slope,classification,actionable,"
        "proposal_method,method_rank,report_selected,loop_diagnostic_status,"
        "loop_artifact_flag,loopless_to_standard_capacity_ratio,loop_diagnostic_reason\n"
        + fseof_rows,
    )
    _write(
        root / "06_amplification" / "fvseof_tidy.csv",
        "target,enforced_product_flux,mean_flux,forced_min_flux,capacity,classification,"
        "park_type,robust,actionable,proposal_method,method_rank,report_selected,"
        "loop_diagnostic_status,loop_artifact_flag,loopless_to_standard_capacity_ratio,"
        "loop_diagnostic_reason\n" + fvseof_rows,
    )

    if many_targets:
        loop_targets = [
            *(f"F{value}" for value in range(1, 11)),
            *(f"V{value}" for value in range(1, 11)),
        ]
    else:
        loop_targets = ["R1", "R2"]
    loop_rows = "".join(
        f"{rank},{target},{'fseof' if target.startswith('F') else 'fvseof'},0,10,10,0,9,9,"
        f"{0.01 if target in {'F1', 'V1'} else 0.9},0.25,"
        f"{'true' if target in {'F1', 'V1'} else 'false'},complete,"
        f"{'loop-dominated diagnostic' if target in {'F1', 'V1'} else ''},2,0.1\n"
        for rank, target in enumerate(loop_targets, start=1)
    )
    _write(
        root / "07_validation" / "amplification_loop_diagnostic.csv",
        "rank,target,source_methods,standard_minimum,standard_maximum,standard_capacity,"
        "loopless_minimum,loopless_maximum,loopless_capacity,"
        "loopless_to_standard_capacity_ratio,capacity_ratio_threshold,loop_artifact_flag,"
        "diagnostic_status,reason,enforced_product_floor,biomass_floor\n" + loop_rows,
    )

    if optional == "complete":
        if exhaustive_validation:
            amplification_targets = (
                [
                    *(f"F{target}" for target in range(1, 11)),
                    *(f"V{target}" for target in range(1, 11)),
                ]
                if many_targets
                else ["R1", "R2"]
            )
            knockout_specs = [
                ("g1", "R1", "g1;g1_alias"),
                ("g2", "R2", "g2;g4"),
                ("g5", "R5", "g5"),
                ("g6", "R6", "g6"),
                ("g7", "R7", "g7"),
                ("g8", "R8", "g8"),
                ("g9", "R9", "g9"),
                ("g10", "R10", "g10"),
                ("g11", "R11", "g11"),
            ]
            response_rows = "".join(
                f"{target},{step},{step + rank / 10},0.3,optimal,{target},EX_product,"
                "wild_type,all_report_selected_candidates\n"
                for rank, target in enumerate(amplification_targets, start=1)
                for step in (0, 4)
            ) + "".join(
                f"{target},{step},{0.2 + step / 2},{0.3 - step / 20},optimal,"
                f"{signature},EX_product,wild_type,all_display_ranked_candidates\n"
                for target, signature, _ in knockout_specs
                for step in (0, 4)
            )
            _write(
                root / "07_validation" / "flux_response_tidy.csv",
                "target,target_flux,response_flux,biomass_flux,status,scan_reaction,"
                "response_reaction,background,candidate_scope\n" + response_rows,
            )
            response_index_header = (
                "target,scan_reaction,response_reaction,background,blocked_reactions,"
                "blocked_reaction_signature,candidate_target_ids,actions,source_methods,"
                "candidate_scope,loop_diagnostic_status,loop_artifact_flag,"
                "loop_diagnostic_eligible,loop_diagnostic_reason,scan_reference_flux,"
                "status,error,reason,"
                "data_file,phases_file,metadata_file\n"
            )
            response_index_rows = "".join(
                f"{target},{target},EX_product,wild_type,,,{target},amplify,fseof,"
                "all_report_selected_candidates,complete,"
                f"{'true' if target in {'F1', 'V1'} else 'false'},"
                f"{'false' if target in {'F1', 'V1'} else 'true'},"
                f"{'loop artifact' if target in {'F1', 'V1'} else ''},1.0,complete,,,"
                f"flux_response__{target}.csv,flux_response_phases__{target}.csv,"
                f"flux_response__{target}.metadata.json\n"
                for target in amplification_targets
            ) + "".join(
                f"{target},{signature},EX_product,wild_type,{signature},{signature},"
                f"{candidate_ids},knockout,moma_l2;room,all_display_ranked_candidates,"
                f",,,,{1e-12 if target == 'g1' else 1.0},complete,,,"
                f"flux_response__{target}.csv,"
                f"flux_response_phases__{target}.csv,flux_response__{target}.metadata.json\n"
                for target, signature, candidate_ids in knockout_specs
            )
            _write(
                root / "07_validation" / "flux_response_index.csv",
                response_index_header + response_index_rows,
            )
            sampling_rows = "".join(
                f"{target},{condition},{reaction},{value},{sample_id}\n"
                for target, _, _ in knockout_specs
                for condition, offset in (("wild_type", 0.0), ("knockout", 0.5))
                for reaction, base in (("EX_product", 0.2), ("BIOMASS", 0.1))
                for sample_id, value in enumerate((base + offset, base + offset + 0.1))
            )
            _write(
                root / "07_validation" / "sampling_tidy.csv",
                "target,condition,reaction_id,flux,sample_id\n" + sampling_rows,
            )
            sampling_index_header = (
                "target_id,blocked_reactions,blocked_reaction_signature,"
                "candidate_target_ids,source_methods,candidate_scope,status,error,reason,"
                "samples_file,statistics_file,comparison_file,metadata_file\n"
            )
            sampling_index_rows = (
                "wild_type,,,wild_type,,,complete,,,random_sampling__wild_type.csv.gz,"
                "random_sampling_statistics__wild_type.csv,,"
                "random_sampling__wild_type.metadata.json\n"
                + "".join(
                    f"{target},{signature},{signature},{candidate_ids},moma_l2;room,"
                    f"all_display_ranked_candidates,complete,,,"
                    f"random_sampling__{target}.csv.gz,"
                    f"random_sampling_statistics__{target}.csv,"
                    f"random_sampling_comparison__{target}.csv,"
                    f"random_sampling__{target}.metadata.json\n"
                    for target, signature, candidate_ids in knockout_specs
                )
            )
            _write(
                root / "07_validation" / "random_sampling_index.csv",
                sampling_index_header + sampling_index_rows,
            )
        else:
            _write(
                root / "07_validation" / "flux_response_tidy.csv",
                "target,target_flux,response_flux,biomass_flux,status,scan_reaction,response_reaction,background\n"
                "R1,0,0.1,0.4,optimal,R1,EX_product,wild_type\n"
                "R1,2,2.1,0.35,optimal,R1,EX_product,wild_type\n"
                "R2,0,0.1,0.4,optimal,R2,EX_product,wild_type\n"
                "R2,3,2.8,0.30,optimal,R2,EX_product,wild_type\n"
                "g2,0,0.2,0.25,optimal,R2,EX_product,gene_knockout\n"
                "g2,4,3.8,0.18,optimal,R2,EX_product,gene_knockout\n",
            )
            _write(
                root / "07_validation" / "sampling_tidy.csv",
                "target,condition,reaction_id,flux,sample_id\n"
                "g2,wild_type,EX_product,0.1,0\ng2,wild_type,EX_product,0.3,1\n"
                "g2,knockout,EX_product,3.5,0\ng2,knockout,EX_product,4.0,1\n"
                "g2,wild_type,BIOMASS,0.4,0\ng2,wild_type,BIOMASS,0.38,1\n"
                "g2,knockout,BIOMASS,0.25,0\ng2,knockout,BIOMASS,0.27,1\n",
            )
        _write(
            root / "07_validation" / "recommendations.csv",
            "target,type,evidence,verdict,proposal_methods,validation_methods,growth_retained,"
            "product_effect,artifact_flag,reason\n"
            "g2,single_gene_knockout,moma_l2;room,support,moma_l2;room,flux_response;sampling,"
            "0.675,3.8,false,model analyses agree\n"
            f"{'F2' if exhaustive_validation and many_targets else 'R1'},amplification,"
            "fseof;fvseof,support,fseof;fvseof,flux_response;loopless_fva,"
            "0.8,2.7,false,forward response supports direction\n"
            '"R1;R2",multi_knockout,optknock;robustknock,coupled,optknock;robustknock,,0.5,'
            "6.5,false,positive guaranteed product\n",
        )
    else:
        _write(
            root / "07_validation" / "recommendations.csv",
            "target,type,evidence,verdict,proposal_methods,validation_methods,growth_retained,"
            "product_effect,artifact_flag,reason\n"
            "R1,amplification,fseof;fvseof,unavailable,fseof;fvseof,,,unknown,false,"
            "forward validation was unavailable\n",
        )

    for relative, content in (
        ("scripts/production_config.json", json.dumps(config, sort_keys=True) + "\n"),
        ("scripts/reproduce.py", "# reproduce fixture\n"),
        ("scripts/render.py", "# render fixture\n"),
        ("scripts/validate.py", "# validate fixture\n"),
    ):
        _write(root / relative, content)

    supplementary: list[dict[str, object]] = []
    artifacts: dict[str, object] = {}
    role_paths = {
        "provenance": "00_provenance.json",
        "summary": "00_summary.json",
        "model": "model/fixture.xml",
        "wild_type_fluxes": "03_reference/wild_type_fluxes.csv",
        "theoretical_yield": "02_yield/theoretical_yield.csv",
        "production_envelope": "02_yield/production_envelope.csv",
        "reference_fluxes": "03_reference/reference_pfba.csv",
        "single_knockout_moma": "04_single_knockout/moma.csv",
        "single_knockout_room": "04_single_knockout/room.csv",
        "single_knockout_consensus": "04_single_knockout/consensus.csv",
        "gene_knockout_mapping": "04_single_knockout/gene_knockout_mapping.csv",
        "optknock": "05_strain_design/optknock.csv",
        "robustknock": "05_strain_design/robustknock.csv",
        "amplification_target_ranking": "06_amplification/fseof.csv",
        "variability_supported_amplification_targets": "06_amplification/fvseof.csv",
        "fseof_tidy": "06_amplification/fseof_tidy.csv",
        "fvseof_tidy": "06_amplification/fvseof_tidy.csv",
        "amplification_loop_diagnostic": "07_validation/amplification_loop_diagnostic.csv",
        "recommendations": "07_validation/recommendations.csv",
        "reproduction_config": "scripts/production_config.json",
        "reproduce_script": "scripts/reproduce.py",
        "render_script": "scripts/render.py",
        "validate_script": "scripts/validate.py",
    }
    if optional == "complete":
        role_paths.update(
            {
                "flux_response_tidy": "07_validation/flux_response_tidy.csv",
                "sampling_tidy": "07_validation/sampling_tidy.csv",
            }
        )
        if exhaustive_validation:
            role_paths.update(
                {
                    "flux_response_validation_index": "07_validation/flux_response_index.csv",
                    "single_knockout_sampling_validation_index": "07_validation/random_sampling_index.csv",
                }
            )
    else:
        artifacts.update(
            {
                "flux_response_tidy": {
                    "status": "skipped",
                    "reason": "No candidate passed the forward-validation gate.",
                },
                "sampling_tidy": {
                    "status": "failed",
                    "reason": "The sampler did not converge under the requested condition.",
                },
            }
        )
    for role, relative in role_paths.items():
        artifacts[role] = _entry(root, relative, role, supplementary)
    fseof_tidy_entry = artifacts["fseof_tidy"]
    assert isinstance(fseof_tidy_entry, dict)
    fseof_tidy_entry["reason"] = "retained 1 flagged candidate(s) as diagnostic-only"

    # Exercise promotion of stable supporting roles from the supplementary inventory.
    for role, relative in (
        ("preflight_checks", "01_preflight/preflight.csv"),
        ("workflow_configuration", "00_config.json"),
    ):
        entry = _entry(root, relative, role, supplementary)
        record = _record(root / relative, root, role)
        record.update(entry)
        supplementary.append(record)

    manifest = {
        "schema_id": "cmm.production-target-discovery",
        "schema_version": 2,
        "report": {
            "language": "en",
            "title": "Fixture production-engineering report",
            "subtitle": "Defined medium; solver fixture",
            "product_label": "Target product",
        },
        "artifacts": artifacts,
        "supplementary_artifacts": supplementary,
    }
    _write(
        root / "00_manifest.json", json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    return root


def _r_is_ready() -> bool:
    executable = shutil.which("Rscript")
    if executable is None:
        return False
    check = subprocess.run(
        [
            executable,
            "--vanilla",
            "-e",
            'quit(status=if(all(vapply(c("jsonlite","ggplot2","ggrepel","patchwork","svglite","ragg"), requireNamespace, logical(1), quietly=TRUE))) 0 else 1)',
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    return check.returncode == 0


def test_validate_schema_v2_integrity_and_pre_render_alias(tmp_path):
    run = _make_run(tmp_path / "run")
    validated = validate_run(run)
    assert validated.manifest["schema_version"] == 2
    assert validated.artifact("single_knockout_moma").name == "moma.csv"
    assert validated.artifact("preflight_checks").name == "preflight.csv"
    report = validate_production_run(run)
    assert report.valid and report.phase == "pre-render"
    assert "have not been rendered" in " ".join(report.warnings)
    assert report.raise_for_errors().root == run.resolve()


def test_legacy_schema_v2_core_columns_validate_and_render_with_fallbacks(tmp_path):
    run = _make_run(tmp_path / "run")
    manifest_path = run / "00_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    legacy_screen = (
        "target_id,kind,status,objective,product_flux,selected,growth_fraction,"
        "blocked_reactions\n"
        "g1,gene,optimal,0.31,3.2,true,0.775,R1\n"
        "g2,gene,optimal,0.27,4.1,true,0.675,R2\n"
        "g4,gene,optimal,0.28,4,false,0.7,R2\n"
        "g5,gene,optimal,0.26,3,false,0.65,R5\n"
        "g6,gene,optimal,0.24,2.5,false,0.6,R6\n"
        "g7,gene,optimal,0.22,2,false,0.55,R7\n"
        "g3,gene,infeasible,,,false,,\n"
    )
    _replace_artifact(run, manifest, "single_knockout_moma", legacy_screen)
    _replace_artifact(run, manifest, "single_knockout_room", legacy_screen)
    _replace_artifact(
        run,
        manifest,
        "gene_knockout_mapping",
        "gene_id,gene_name,inert,blocked_reaction,reaction_name,gpr\n"
        "g1,Gene 1,false,R1,Reaction 1,g1\n"
        "g2,Gene 2,false,R2,Reaction 2,g2 and g4\n"
        "g4,Gene 4,false,R2,Reaction 2,g2 and g4\n",
    )
    _replace_artifact(
        run,
        manifest,
        "fseof_tidy",
        "target,enforced_product_flux,reaction_flux\n"
        "R1,0,1\nR1,12,5\nR2,0,2\nR2,12,4\n",
    )
    _replace_artifact(
        run,
        manifest,
        "fvseof_tidy",
        "target,enforced_product_flux,mean_flux,forced_min_flux\n"
        "R1,0,1,0.5\nR1,12,5,4.5\nR2,0,2,1\nR2,12,4,3.5\n",
    )
    artifacts = manifest["artifacts"]
    assert isinstance(artifacts, dict)
    del artifacts["amplification_target_ranking"]
    del artifacts["variability_supported_amplification_targets"]
    manifest_path.write_text(json.dumps(manifest))

    validated = validate_run(run)
    assert validated.manifest["schema_version"] == 2
    assert "amplification_target_ranking" not in validated.artifacts
    if _r_is_ready():
        figures = render_publication_figures(run)
        statuses = {item["id"]: item["status"] for item in figures.figures}
        assert statuses["fig02_single_knockout"] == "rendered"
        assert statuses["fig04_amplification"] == "rendered"
        knockout_svg = (run / "figures/fig02_single_knockout.svg").read_text()
        assert knockout_svg.count(">D1 ") == 2
        report = build_publication_report(run).report_html.read_text()
        assert "FSEOF independent top ten" in report
        assert "R1" in report


def test_validator_aggregates_schema_integrity_and_metadata_errors(tmp_path):
    run = _make_run(tmp_path / "run")
    manifest_path = run / "00_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    del manifest["artifacts"]["single_knockout_room"]
    manifest["artifacts"]["single_knockout_moma"]["sha256"] = "0" * 64
    del manifest["artifacts"]["production_envelope"]["metadata_path"]
    manifest_path.write_text(json.dumps(manifest))

    with pytest.raises(RunValidationError) as raised:
        validate_run(run)
    message = str(raised.value)
    assert "single_knockout_room" in message
    assert "SHA-256 does not match" in message
    assert "must declare metadata_path" in message
    result = validate_production_run(run)
    assert not result.valid and len(result.issues) >= 3


def test_amplification_negative_result_uses_independent_method_support(tmp_path):
    run = _make_run(tmp_path / "run", many_targets=True)
    manifest_path = run / "00_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    _replace_artifact(
        run,
        manifest,
        "recommendations",
        "target,type,evidence,verdict,proposal_methods,validation_methods,growth_retained,"
        "product_effect,artifact_flag,reason\n",
    )
    manifest_path.write_text(json.dumps(manifest))
    validated = validate_run(run)
    statement = _amplification_support_statement(validated, ())
    assert "No independently proposed FSEOF or FVSEOF target" in statement
    assert "method-specific forward-validation support rule" in statement
    assert "shared target" not in statement


@pytest.mark.parametrize(
    ("mutation", "match"),
    [
        (lambda manifest: manifest.__setitem__("schema_version", 2.0), "integer 2"),
        (
            lambda manifest: manifest["artifacts"]["model"].__setitem__(
                "path", "../outside.xml"
            ),
            "stay inside",
        ),
        (
            lambda manifest: manifest["artifacts"]["model"].__setitem__(
                "size_bytes", -1
            ),
            "non-negative integer",
        ),
    ],
)
def test_validator_rejects_manifest_boundary_errors(tmp_path, mutation, match):
    run = _make_run(tmp_path / "run")
    manifest_path = run / "00_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    mutation(manifest)
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(RunValidationError, match=match):
        validate_run(run)


def test_validator_rejects_scientifically_inconsistent_design(tmp_path):
    run = _make_run(tmp_path / "run")
    path = run / "05_strain_design" / "optknock.csv"
    path.write_text(
        "knockouts,growth,max_product,guaranteed_product,growth_coupled\nR1,0.2,2,3,true\n"
    )
    manifest_path = run / "00_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    entry = manifest["artifacts"]["optknock"]
    entry["sha256"] = _digest(path)
    entry["size_bytes"] = path.stat().st_size
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(
        RunValidationError, match="guaranteed_product above max_product"
    ):
        validate_run(run)


def test_validator_aggregates_numeric_and_cross_method_scientific_failures(tmp_path):
    run = _make_run(tmp_path / "run")
    manifest_path = run / "00_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    _replace_artifact(
        run,
        manifest,
        "summary",
        json.dumps({"status": "partial", "model_sha256": "b" * 64}) + "\n",
    )
    _replace_artifact(
        run,
        manifest,
        "wild_type_fluxes",
        "reaction_id,flux\n,NaN\nR1,1\nR1,2\n",
    )
    _replace_artifact(
        run,
        manifest,
        "theoretical_yield",
        "status,molar_yield,product_flux,substrate_uptake,carbon_imbalance\n"
        "infeasible,-1,0,NaN,false\n",
    )
    _replace_artifact(
        run,
        manifest,
        "production_envelope",
        "product_flux,growth_min,growth_max\n2,-1,-2\n1,0,0\n",
    )
    screen_header = (
        "target_id,kind,status,objective_value,distance,distance_kind,objective,"
        "n_reactions,product_flux,selected\n"
    )
    _replace_artifact(
        run,
        manifest,
        "single_knockout_moma",
        screen_header + ",gene,optimal,1,1,l1,-0.1,1,-2,true\n",
    )
    _replace_artifact(
        run,
        manifest,
        "single_knockout_room",
        screen_header + "different,gene,optimal,1,1,l1,0.1,1,2,true\n",
    )
    design_header = "knockouts,growth,max_product,guaranteed_product,growth_coupled\n"
    _replace_artifact(
        run,
        manifest,
        "optknock",
        design_header + "R1,-1,-2,-3,maybe\nR2,0.1,2,3,true\n",
    )
    _replace_artifact(
        run,
        manifest,
        "fseof_tidy",
        "target,enforced_product_flux,reaction_flux,classification\nR1,NaN,1,amplification\n",
    )
    _replace_artifact(
        run,
        manifest,
        "fvseof_tidy",
        "target,enforced_product_flux,mean_flux,forced_min_flux,classification,robust\n"
        "R1,0,1,,amplification,true\n",
    )
    _replace_artifact(
        run,
        manifest,
        "flux_response_tidy",
        "target,target_flux,response_flux,biomass_flux,status,scan_reaction,response_reaction,background\n"
        ",NaN,NaN,,optimal,,,other\n",
    )
    _replace_artifact(
        run,
        manifest,
        "sampling_tidy",
        "target,condition,reaction_id,flux,sample_id\n"
        "g2,wild_type,EX_product,NaN,0\ng3,other,EX_product,1,0\n",
    )
    loop_header = (
        "rank,target,source_methods,standard_minimum,standard_maximum,standard_capacity,"
        "loopless_minimum,loopless_maximum,loopless_capacity,"
        "loopless_to_standard_capacity_ratio,capacity_ratio_threshold,loop_artifact_flag,"
        "diagnostic_status,reason,enforced_product_floor,biomass_floor\n"
    )
    _replace_artifact(
        run,
        manifest,
        "amplification_loop_diagnostic",
        loop_header + "1,R1,fseof,,,,,,,,,maybe,complete,,,,\n",
    )
    _replace_artifact(
        run,
        manifest,
        "recommendations",
        "target,type,evidence,verdict,proposal_methods,validation_methods,growth_retained,"
        "product_effect,artifact_flag,reason\n"
        ",invalid,x,support,x,x,true,increase,false,bad target\n"
        "X,amplification,x,unexpected,x,x,true,increase,false,bad verdict\n"
        "Z,single_gene_knockout,x,support,x,x,true,increase,false,no validation\n"
        "A,amplification,x,support,x,x,true,increase,false,no validation\n",
    )
    manifest_path.write_text(json.dumps(manifest))

    with pytest.raises(RunValidationError) as raised:
        validate_run(run)
    message = str(raised.value)
    for expected in (
        "summary status must be 'complete'",
        "summary model_sha256 does not match",
        "repeats reaction_id",
        "theoretical_yield status",
        "growth_min greater than growth_max",
        "product_flux must be non-decreasing",
        "MOMA and ROOM must cover the same knockout universe",
        "invalid growth_coupled flag",
        "fewer than two scan points",
        "invalid background",
        "matched wild_type and knockout",
        "invalid loop_artifact_flag",
        "invalid verdict",
        "no paired sampling rows",
    ):
        assert expected in message


@pytest.fixture(scope="module")
def rendered_run(tmp_path_factory):
    if not _r_is_ready():
        pytest.skip("Rscript publication packages are not installed")
    run = _make_run(tmp_path_factory.mktemp("publication") / "run")
    bundle = render_production_report(run)
    return run, bundle


def test_r_renderer_writes_validated_vector_and_300_dpi_artwork(rendered_run):
    run, bundle = rendered_run
    rendered = {
        item["id"]: item
        for item in bundle.figures.figures
        if item["status"] == "rendered"
    }
    assert set(rendered) == {
        "fig01_yield_envelope",
        "fig02_single_knockout",
        "fig03_strain_design",
        "fig04_amplification",
        "fig05_flux_response",
        "fig06_sampling_shift",
    }
    check = validate_production_run(run)
    assert check.valid and check.phase == "post-render"
    artwork = check.checks["artwork"]
    assert artwork["status"] == "pass"
    assert len(artwork["files"]) == 18
    pngs = [item for item in artwork["files"] if item["format"] == "png"]
    assert all(item["dpi_x"] == pytest.approx(300, abs=1) for item in pngs)
    assert bundle.report.report_validation == run / "report_validation.json"
    assert bundle.figures.renderer["script_sha256"] == _digest(renderer_script_path())

    figure_two = rendered["fig02_single_knockout"]
    assert figure_two["sources"] == [
        "04_single_knockout/moma.csv",
        "04_single_knockout/room.csv",
        "00_summary.json",
        "scripts/production_config.json",
        "04_single_knockout/consensus.csv",
        "07_validation/recommendations.csv",
    ]
    assert "Every D1-D5 row is a forward-validation candidate" in figure_two["caption"]
    assert "final single-gene support" in figure_two["caption"]
    knockout_svg = (run / "figures/fig02_single_knockout.svg").read_text()
    assert "Validation candidate (D1-D5)" in knockout_svg
    assert "Supported recommendation" in knockout_svg
    assert "Wild-type reference" in knockout_svg
    for rank in range(1, 6):
        assert knockout_svg.count(f">D{rank} ") == 2
    assert knockout_svg.count(">D3 g5</text>") == 2
    assert knockout_svg.count(">D5 g7</text>") == 2
    assert "not a recommendation" in figure_two["caption"]

    figure_four = rendered["fig04_amplification"]
    assert "retained in flux-response validation" in figure_four["caption"]
    assert (
        "excluded from recommendation and forward validation"
        not in figure_four["caption"]
    )
    amplification_svg = (run / "figures/fig04_amplification.svg").read_text()
    direct_label_styles = re.findall(
        r"<text[^>]+style='([^']+)'[^>]*>D\d+ [^<]+</text>",
        amplification_svg,
    )
    assert direct_label_styles
    for style in direct_label_styles:
        font_size = re.search(r"font-size:\s*([0-9.]+)px", style)
        fill = re.search(r"fill:\s*(#[0-9A-Fa-f]{6})", style)
        assert font_size is not None and float(font_size.group(1)) >= 5.0
        assert fill is not None
        channels = [
            int(fill.group(1)[offset : offset + 2], 16) / 255 for offset in (1, 3, 5)
        ]
        linear = [
            value / 12.92 if value <= 0.04045 else ((value + 0.055) / 1.055) ** 2.4
            for value in channels
        ]
        luminance = 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]
        white_contrast = 1.05 / (luminance + 0.05)
        assert white_contrast >= 4.5

    figure_three = rendered["fig03_strain_design"]
    assert "guaranteed endpoint is drawn on top" in figure_three["caption"]

    figure_five = rendered["fig05_flux_response"]
    assert "Standard target-to-product flux-response curves" in figure_five["caption"]
    assert "enforced candidate-reaction flux (target_flux)" in figure_five["caption"]
    assert "optimized target-product flux (response_flux)" in figure_five["caption"]
    assert "not a plot axis" in figure_five["caption"]
    assert (
        "Legacy unscoped panels retain their recorded model background"
        in figure_five["caption"]
    )
    assert (
        "enforced-candidate-reaction-flux versus target-product-flux"
        in figure_five["alt"]
    )
    assert "target(s)" not in figure_five["caption"]
    assert "background(s)" not in figure_five["caption"]
    response_svg = (run / "figures/fig05_flux_response.svg").read_text()
    assert "Amplification candidates (wild type)" in response_svg
    assert (
        "Knockout-derived candidate reaction scans (recorded backgrounds)"
        in response_svg
    )
    assert ">g2 (R2)</text>" in response_svg
    assert response_svg.count(">Enforced candidate-reaction flux</text>") == 2
    assert response_svg.count(">Target-product flux</text>") == 2
    assert ">Growth rate</text>" not in response_svg
    assert "EX_product" not in response_svg
    assert "BIOMASS" not in response_svg
    renderer_source = renderer_script_path().read_text()
    assert (
        renderer_source.count(
            "ggplot2::aes(x = target_flux, y = response_flux, group = target)"
        )
        == 2
    )
    assert renderer_source.count("limits = shared_response_limits") == 2
    assert renderer_source.count('scales = "free_x"') == 2
    assert renderer_source.count("min(4L, length(") >= 2
    assert renderer_source.count("scale_x_continuous(breaks = facet_flux_breaks)") == 2
    assert "pretty(limits, n = 3)" in renderer_source
    assert renderer_source.count('panel.spacing.x = grid::unit(4, "mm")') == 2
    assert "at most four facets per row" in figure_five["caption"]
    assert "aes(x = biomass_flux" not in renderer_source

    figure_six = rendered["fig06_sampling_shift"]
    assert "intentionally omits other reactions" in figure_six["caption"]
    sampling_svg = (run / "figures/fig06_sampling_shift.svg").read_text()
    assert "Product exchange" in sampling_svg
    assert "Biomass reaction" in sampling_svg
    assert "EX_product" not in sampling_svg
    assert "BIOMASS" not in sampling_svg


def test_exhaustive_candidate_validation_is_complete_and_publication_readable(tmp_path):
    if not _r_is_ready():
        pytest.skip("Rscript publication packages are not installed")
    run = _make_run(tmp_path / "run", many_targets=True, exhaustive_validation=True)
    validated = validate_run(run)
    response_index = publication._rows(validated, "flux_response_validation_index")
    sampling_index = publication._rows(
        validated, "single_knockout_sampling_validation_index"
    )
    assert len(response_index) == 29
    assert len([row for row in sampling_index if row["target_id"] != "wild_type"]) == 9

    bundle = render_production_report(run)
    rendered = {
        item["id"]: item
        for item in bundle.figures.figures
        if item["status"] == "rendered"
    }
    figure_five = rendered["fig05_flux_response"]
    figure_six = rendered["fig06_sampling_shift"]
    assert figure_five["height_mm"] == 230
    assert figure_six["height_mm"] == 174
    assert figure_five["sources"] == [
        "07_validation/flux_response_tidy.csv",
        "07_validation/flux_response_index.csv",
    ]
    assert figure_six["sources"] == [
        "07_validation/sampling_tidy.csv",
        "00_summary.json",
        "07_validation/random_sampling_index.csv",
    ]

    response_svg = (run / "figures/fig05_flux_response.svg").read_text()
    for target in [
        *(f"F{value}" for value in range(1, 11)),
        *(f"V{value}" for value in range(1, 11)),
    ]:
        assert f">{target}</text>" in response_svg
    for target in ("g1", "g2", "g5", "g6", "g7", "g8", "g9", "g10", "g11"):
        reaction = target.replace("g", "R")
        assert f">{target} ({reaction})</text>" in response_svg
    sampling_svg = (run / "figures/fig06_sampling_shift.svg").read_text()
    for target in ("g1", "g2", "g5", "g6", "g7", "g8", "g9", "g10", "g11"):
        assert sampling_svg.count(f">{target}</text>") == 2

    linked = bundle.report.report_html.read_text()
    assert "20/20" in linked
    assert linked.count("9/9") >= 2
    assert "29" in figure_five["caption"]
    assert "wild-type pre-deletion single-reaction titrations" in figure_five["caption"]
    assert "representative gene with its scanned reaction" in figure_five["caption"]
    assert "1 zero-reference knockout-derived scan (g1)" in figure_five["caption"]
    assert "cannot causally support deletion" in figure_five["caption"]
    assert "all 9 completed" in figure_six["caption"]
    assert "g1;g1_alias" in linked
    assert "share one model intervention" in linked
    assert "Candidate coverage is separate from recommendation status" in linked
    assert "Zero-reference knockout candidates" in linked
    assert "full feasible candidate-reaction domain" in linked
    assert "Reference candidate flux (mmol gDW⁻¹ h⁻¹)" in linked
    assert "Legacy schema-v2 coverage" not in linked
    check = validate_production_run(run)
    assert check.valid and check.phase == "post-render"


@pytest.mark.parametrize(
    ("role", "relative_path", "target_column"),
    [
        (
            "flux_response_validation_index",
            "07_validation/flux_response_index.csv",
            "target",
        ),
        (
            "single_knockout_sampling_validation_index",
            "07_validation/random_sampling_index.csv",
            "target_id",
        ),
    ],
)
def test_exhaustive_indexes_require_every_signature_equivalent_gene_alias(
    tmp_path, role, relative_path, target_column
):
    run = _make_run(tmp_path / "run", many_targets=True, exhaustive_validation=True)
    manifest_path = run / "00_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    path = run / relative_path
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
        fields = list(rows[0])
    target_row = next(row for row in rows if row[target_column] == "g1")
    assert target_row["candidate_target_ids"] == "g1;g1_alias"
    target_row["candidate_target_ids"] = "g1"
    content = (
        ",".join(fields)
        + "\n"
        + "".join(",".join(row[field] for field in fields) + "\n" for row in rows)
    )
    _replace_artifact(run, manifest, role, content)
    manifest_path.write_text(json.dumps(manifest))

    result = validate_production_run(run)
    assert not result.valid
    assert "must list every signature-equivalent gene" in " ".join(result.issues)
    assert "g1_alias" in " ".join(result.issues)


def test_current_flux_response_scope_requires_wild_type_background(tmp_path):
    run = _make_run(tmp_path / "run", many_targets=True, exhaustive_validation=True)
    manifest_path = run / "00_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    for role in ("flux_response_validation_index", "flux_response_tidy"):
        _update_csv_artifact_rows(
            run,
            manifest,
            role,
            target_column="target",
            target="g1",
            updates={"background": "gene_knockout"},
        )
    manifest_path.write_text(json.dumps(manifest))

    result = validate_production_run(run)
    assert not result.valid
    assert "current candidate_scope" in " ".join(result.issues)
    assert "background is not 'wild_type'" in " ".join(result.issues)


def test_current_flux_response_scope_requires_configured_product_response(tmp_path):
    run = _make_run(tmp_path / "run", many_targets=True, exhaustive_validation=True)
    manifest_path = run / "00_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    for role in ("flux_response_validation_index", "flux_response_tidy"):
        _update_csv_artifact_rows(
            run,
            manifest,
            role,
            target_column="target",
            target="g1",
            updates={"response_reaction": "BIOMASS"},
        )
    manifest_path.write_text(json.dumps(manifest))

    result = validate_production_run(run)
    assert not result.valid
    assert "must use configured product 'EX_product'" in " ".join(result.issues)


def test_current_flux_response_scope_requires_finite_secondary_biomass(tmp_path):
    run = _make_run(tmp_path / "run", many_targets=True, exhaustive_validation=True)
    manifest_path = run / "00_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    _update_csv_artifact_rows(
        run,
        manifest,
        "flux_response_tidy",
        target_column="target",
        target="g1",
        updates={"biomass_flux": ""},
    )
    manifest_path.write_text(json.dumps(manifest))

    result = validate_production_run(run)
    assert not result.valid
    assert "has no finite 'biomass_flux' value" in " ".join(result.issues)


def test_current_flux_response_index_cannot_be_downgraded_to_unscoped_tidy(tmp_path):
    run = _make_run(tmp_path / "run", many_targets=True, exhaustive_validation=True)
    manifest_path = run / "00_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    _update_csv_artifact_rows(
        run,
        manifest,
        "flux_response_tidy",
        target_column="target",
        target="g1",
        updates={
            "candidate_scope": "",
            "scan_reaction": "WRONG_SCAN",
            "response_reaction": "WRONG_RESPONSE",
            "background": "gene_knockout",
            "biomass_flux": "",
        },
    )
    manifest_path.write_text(json.dumps(manifest))

    result = validate_production_run(run)
    assert not result.valid
    assert "has no tidy rows with matching target 'g1'" in " ".join(result.issues)
    assert "candidate_scope 'all_display_ranked_candidates'" in " ".join(result.issues)


def test_exhaustive_indexes_report_failed_and_skipped_candidates(tmp_path):
    run = _make_run(tmp_path / "run", many_targets=True, exhaustive_validation=True)
    manifest_path = run / "00_manifest.json"
    manifest = json.loads(manifest_path.read_text())

    response_index_path = run / "07_validation" / "flux_response_index.csv"
    with response_index_path.open(encoding="utf-8", newline="") as handle:
        response_rows = list(csv.DictReader(handle))
        response_fields = list(response_rows[0])
    for row in response_rows:
        if row["target"] == "V10":
            row["status"] = "failed"
            row["error"] = "SolverError: response solve timed out"
        elif row["target"] == "g11":
            row["status"] = "skipped"
            row["reason"] = (
                "multi-reaction blocked signature has no single candidate-reaction scan"
            )
            row["data_file"] = ""
            row["phases_file"] = ""
            row["metadata_file"] = ""
    response_content = (
        ",".join(response_fields)
        + "\n"
        + "".join(
            ",".join(row[field] for field in response_fields) + "\n"
            for row in response_rows
        )
    )
    _replace_artifact(run, manifest, "flux_response_validation_index", response_content)
    response_tidy_path = run / "07_validation" / "flux_response_tidy.csv"
    response_lines = response_tidy_path.read_text().splitlines()
    _replace_artifact(
        run,
        manifest,
        "flux_response_tidy",
        "\n".join(
            line
            for index, line in enumerate(response_lines)
            if index == 0 or not (line.startswith("V10,") or line.startswith("g11,"))
        )
        + "\n",
    )

    sampling_index_path = run / "07_validation" / "random_sampling_index.csv"
    with sampling_index_path.open(encoding="utf-8", newline="") as handle:
        sampling_rows = list(csv.DictReader(handle))
        sampling_fields = list(sampling_rows[0])
    for row in sampling_rows:
        if row["target_id"] == "g11":
            row["status"] = "skipped"
            row["reason"] = "shared wild-type sampler was unavailable"
    sampling_content = (
        ",".join(sampling_fields)
        + "\n"
        + "".join(
            ",".join(row[field] for field in sampling_fields) + "\n"
            for row in sampling_rows
        )
    )
    _replace_artifact(
        run,
        manifest,
        "single_knockout_sampling_validation_index",
        sampling_content,
    )
    sampling_tidy_path = run / "07_validation" / "sampling_tidy.csv"
    sampling_lines = sampling_tidy_path.read_text().splitlines()
    _replace_artifact(
        run,
        manifest,
        "sampling_tidy",
        "\n".join(
            line
            for index, line in enumerate(sampling_lines)
            if index == 0 or not line.startswith("g11,")
        )
        + "\n",
    )
    summary_payload = json.loads((run / "00_summary.json").read_text())
    summary_payload["validation_coverage"].update(
        {
            "flux_response_attempted": 28,
            "flux_response_completed": 27,
            "flux_response_failed": 1,
            "sampling_attempted": 9,
            "sampling_completed": 9,
            "sampling_skipped": 1,
        }
    )
    _replace_artifact(
        run,
        manifest,
        "summary",
        json.dumps(summary_payload, sort_keys=True) + "\n",
    )
    manifest_path.write_text(json.dumps(manifest))

    validated = validate_run(run)
    summary = publication._validation_execution_summary(validated)
    assert "20/20" in summary
    assert summary.count("9/9") >= 2
    assert "SolverError: response solve timed out" in summary
    assert (
        "multi-reaction blocked signature has no single candidate-reaction scan"
        in summary
    )
    assert "shared wild-type sampler was unavailable" in summary
    assert "Legacy schema-v2 coverage" not in summary


def test_html_is_structured_deterministic_linked_and_standalone(rendered_run):
    run, _ = rendered_run
    first = build_publication_report(run)
    linked_once = first.report_html.read_bytes()
    standalone_once = first.report_standalone_html.read_bytes()
    second = build_publication_report(run)
    assert second.report_html.read_bytes() == linked_once
    assert second.report_standalone_html.read_bytes() == standalone_once

    linked = linked_once.decode()
    standalone = standalone_once.decode()
    for heading in (
        "1. Summary",
        "2. Setup",
        "3. Data and methods",
        "4. Results",
        "5. Recommended targets and strain proposal",
        "6. Limitations",
        "7. References",
        "8. Provenance",
    ):
        assert heading in linked
    assert 'src="figures/fig01_yield_envelope.png"' in linked
    assert 'style="width:89mm;max-width:100%"' in linked
    assert 'style="width:180mm;max-width:100%"' in linked
    assert 'href="02_yield/theoretical_yield.csv"' in linked
    assert "data:image/png;base64," in standalone
    assert 'href="02_yield/' not in standalone
    assert str(run.resolve()) not in linked
    assert str(run.resolve()) not in standalone
    assert "Run directory: <code>.</code>" in linked
    assert "fixture condition warning" in linked
    assert linked.index("fixture condition warning") < linked.index("2. Setup")
    assert "summary oxygen-bound warning" in linked
    assert linked.index("summary oxygen-bound warning") < linked.index("2. Setup")
    summary_section = linked[: linked.index("2. Setup")]
    assert summary_section.count("medium was used as loaded") == 1
    assert (
        "model bounds were retained under the model-as-loaded medium"
        not in summary_section
    )
    assert ".preflight-table th:nth-child(2)" in linked
    assert "white-space:nowrap" in linked
    assert "Gene-to-reaction interpretation" in linked
    assert "Loopless capacity diagnostic" in linked
    assert (
        "Optimized response" not in linked
    )  # axis text lives in vector/raster artwork
    assert "confidence score" in linked
    assert "https://doi.org/10.1093/bioinformatics/btp704" in linked
    assert "https://doi.org/10.1093/bioinformatics/btac632" in linked
    assert "https://doi.org/10.1093/bioinformatics/btw555" in linked
    assert "https://doi.org/10.1074/jbc.R800048200" in linked
    assert "single-gene knockout g2" in linked
    assert "67.5% WT growth retained" in linked
    assert "product-flux shift +3.8 mmol gDW⁻¹ h⁻¹" in linked
    assert "amplification target R1" in linked
    assert "1 growth-coupled reaction-level multi-knockout design" in linked
    assert "R1;R2" in linked
    assert "50% WT growth retained" in linked
    assert "WT growth retained (%)" in linked
    assert "67.5%" in linked
    assert "Conservative sampling Δ product flux" in linked
    assert "Response Δ product flux" in linked
    assert "Guaranteed product flux" in linked
    assert "sampling figure is restricted" in linked
    assert "not independent causal evidence" in linked
    assert "standard candidate-reaction-to-product flux-response definition" in linked
    assert "target_flux</code>) is on the x-axis" in linked
    assert "response_flux</code>) is on the y-axis" in linked
    assert "Biomass flux is a secondary value" in linked
    assert "Legacy schema-v2 rows retain their" in linked
    assert "recorded model background" in linked
    assert "multi-reaction knockout signature remains" in linked
    assert "relabelled as product responses" in linked
    assert "Signature-equivalent rows are model-equivalent deletions" in linked
    assert "Signature-equivalent model deletion" in linked
    assert "g_cross" not in linked
    assert "g4" in linked
    assert "Blocked reactions" in linked
    assert "Reaction equations" in linked
    assert 'class="table-scroll"' in linked
    assert 'class="ko-screen-table"' in linked
    assert "g2 (Gene 2)" in linked
    assert "R2 — Reaction 2" in linked
    assert "R2: B --&gt; C" in linked
    assert "Display rank is not a recommendation" in linked
    assert linked.count("all 5 canonical candidates") == 2
    assert "Beneficial screen candidate (forward-tested)" in linked
    assert "Display candidate (forward-tested)" in linked
    assert "Validated beneficial candidate" not in linked
    assert "Display only" not in linked
    assert "Method shortlist" not in linked
    assert "Every D1–D5 display-ranked signature is a" in linked
    assert "forward-validation shortlist" not in linked
    assert "FSEOF independent top ten" in linked
    assert "FVSEOF independent top ten" in linked
    assert "intersection is not required" in linked
    assert "retained in flux-response validation" in linked
    assert "not eligible for recommendation or forward validation" not in linked
    assert "ineligible for recommendation or forward validation" not in linked
    assert "excluded from recommendation and forward validation" not in linked
    assert "retained 1 flagged candidate as diagnostic-only" in linked
    figure_manifest = (run / "figures/figure_manifest.json").read_text()
    assert "(s)" not in linked
    assert "(s)" not in standalone
    assert "(s)" not in figure_manifest
    assert (
        "Medium: Applied</td><td>{&quot;EX_o2&quot;: 0.0, &quot;EX_substrate&quot;: 10.0}</td>"
        in linked
    )


def test_optional_panels_are_explicit_without_blocking_source_validation(tmp_path):
    if not _r_is_ready():
        pytest.skip("Rscript publication packages are not installed")
    run = _make_run(tmp_path / "run", optional="skipped")
    manifest = render_publication_figures(run)
    statuses = {item["id"]: item["status"] for item in manifest.figures}
    assert statuses["fig01_yield_envelope"] == "rendered"
    assert statuses["fig05_flux_response"] == "skipped"
    assert statuses["fig06_sampling_shift"] == "skipped"
    report = build_publication_report(run)
    document = report.report_html.read_text()
    assert "Validation panel unavailable" in document
    assert "sampler did not converge" in document
    assert validate_production_run(run).valid


def test_r_smoke_handles_empty_designs_and_independent_top_ten_amplification_targets(
    tmp_path,
):
    if not _r_is_ready():
        pytest.skip("Rscript publication packages are not installed")
    run = _make_run(tmp_path / "run", empty_design=True, many_targets=True)
    manifest = render_publication_figures(run)
    statuses = {item["id"]: item["status"] for item in manifest.figures}
    assert statuses["fig03_strain_design"] == "rendered"
    assert statuses["fig04_amplification"] == "rendered"
    assert (run / "figures/fig03_strain_design.svg").stat().st_size > 1000
    figure_four = run / "figures/fig04_amplification.svg"
    svg = figure_four.read_text()
    assert "Target-reaction flux" in svg
    assert "FSEOF: independent top ten" in svg
    assert "FVSEOF: independent top ten" in svg
    assert svg.count("Loop diagnostic only") == 2
    assert svg.count("Eligible top-ranked targets") == 2
    assert "D1 F1 [loop]" in svg
    assert "D1 V1 [loop]" in svg
    for target in [
        *(f"F{value}" for value in range(1, 11)),
        *(f"V{value}" for value in range(1, 11)),
    ]:
        assert target in svg


def test_renderer_declares_compatible_version_floors_and_portable_utf8_environment():
    script = renderer_script_path().read_text(encoding="utf-8")
    for required in (
        'minimum_r_version <- "4.3.2"',
        'jsonlite = "1.8.8"',
        'ggplot2 = "3.5.2"',
        'ggrepel = "0.9.6"',
        'patchwork = "1.2.0"',
        'svglite = "2.2.1"',
        'ragg = "1.2.7"',
        "actual %s; required >= %s",
        "Restore the checked-in renv.lock",
    ):
        assert required in script

    windows = _renderer_environment("a" * 64, platform_name="nt")
    posix = _renderer_environment("a" * 64, platform_name="posix")
    assert windows["LANG"] == "C.UTF-8"
    assert "LC_ALL" not in windows
    assert posix["LC_ALL"] == "C.UTF-8"
    assert windows["CMM_RENDERER_SHA256"] == "a" * 64
    assert _decode_renderer_stream("target → product".encode()) == "target → product"
    assert "�" in _decode_renderer_stream(b"bad:\xff")


def test_renderer_rejects_warning_even_with_zero_exit(tmp_path):
    if not _r_is_ready():
        pytest.skip("Rscript publication packages are not installed")
    run = _make_run(tmp_path / "run")
    script = tmp_path / "warning_renderer.R"
    _write(script, 'warning("fixture plotting warning")\n')
    with pytest.raises(FigureRenderError, match="emitted a warning"):
        render_publication_figures(run, renderer=script)


def test_renderer_decodes_subprocess_bytes_explicitly(tmp_path, monkeypatch):
    run = _make_run(tmp_path / "run")
    captured: dict[str, object] = {}

    def fake_run(command, **kwargs):
        captured["command"] = command
        captured.update(kwargs)
        return subprocess.CompletedProcess(
            command,
            returncode=1,
            stdout="renderer µ".encode(),
            stderr=b"invalid byte: \xff",
        )

    monkeypatch.setattr(publication.shutil, "which", lambda _: "/fake/Rscript")
    monkeypatch.setattr(publication.subprocess, "run", fake_run)
    with pytest.raises(FigureRenderError) as raised:
        render_publication_figures(run)
    assert "renderer µ" in str(raised.value)
    assert "invalid byte: �" in str(raised.value)
    assert "text" not in captured
    environment = captured["env"]
    assert isinstance(environment, dict)
    assert environment["CMM_RENDERER_SHA256"] == _digest(renderer_script_path())


def test_post_render_validator_detects_raster_embedded_svg_and_stale_record(
    rendered_run, tmp_path
):
    source, _ = rendered_run
    run = tmp_path / "tampered"
    shutil.copytree(source, run)
    svg = run / "figures/fig01_yield_envelope.svg"
    document = svg.read_text()
    svg.write_text(
        document.replace("</svg>", '<image href="data:image/png;base64,AA=="/></svg>')
    )

    result = validate_production_run(run)
    assert not result.valid
    assert "embeds a raster image" in " ".join(result.issues)
    assert "stale" in " ".join(result.issues)


def test_post_render_validator_detects_horizontal_svg_text_clipping(
    rendered_run, tmp_path
):
    source, _ = rendered_run
    run = tmp_path / "tampered"
    shutil.copytree(source, run)
    svg = run / "figures/fig04_amplification.svg"
    document = svg.read_text()
    svg.write_text(
        document.replace(
            "</svg>",
            "<text x='520' y='10' textLength='20px'>clipped legend</text></svg>",
        )
    )

    result = validate_production_run(run)
    assert not result.valid
    assert "horizontally clipped text: clipped legend" in " ".join(result.issues)
    assert "stale" in " ".join(result.issues)


def test_post_render_validator_detects_broken_link_and_standalone_relative_href(
    rendered_run, tmp_path
):
    source, _ = rendered_run
    run = tmp_path / "tampered"
    shutil.copytree(source, run)
    linked = run / "report.html"
    linked.write_text(
        linked.read_text().replace("02_yield/theoretical_yield.csv", "missing.csv")
    )
    standalone = run / "report_standalone.html"
    standalone.write_text(
        standalone.read_text().replace(
            "</main>", '<a href="figures/local.csv">bad</a></main>'
        )
    )

    result = validate_production_run(run)
    assert not result.valid
    assert "href does not exist" in " ".join(result.issues)
    assert "relative href" in " ".join(result.issues)


def test_post_render_validator_rejects_absolute_local_path_leak(rendered_run, tmp_path):
    source, _ = rendered_run
    run = tmp_path / "tampered"
    shutil.copytree(source, run)
    standalone = run / "report_standalone.html"
    standalone.write_text(
        standalone.read_text().replace("</main>", f"<p>{run.resolve()}</p></main>")
    )

    result = validate_production_run(run)
    assert not result.valid
    assert "absolute local run directory" in " ".join(result.issues)


def test_build_report_fails_strictly_when_vector_output_is_missing(
    rendered_run, tmp_path
):
    source, _ = rendered_run
    run = tmp_path / "tampered"
    shutil.copytree(source, run)
    (run / "figures/fig02_single_knockout.pdf").unlink()
    with pytest.raises(FigureRenderError, match="missing non-empty pdf"):
        build_publication_report(run)
