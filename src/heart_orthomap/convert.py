"""Sparse, no-duplication count projections; raw inputs are never modified."""
from copy import deepcopy
from dataclasses import dataclass
import hashlib
import re

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .mapping import MappingTable, _canonical


@dataclass
class ConversionResult:
    adata: ad.AnnData
    decisions: pd.DataFrame
    report: dict


def raw_counts(matrix):
    if sparse.issparse(matrix):
        if not getattr(matrix, "has_canonical_format", True):
            raise ValueError("Canonical sparse counts required; resolve duplicate sparse entries first")
        values = matrix.data
    else:
        matrix = np.asarray(matrix)
        if matrix.ndim != 2:
            raise ValueError("Expected a cells-by-genes matrix")
        values = matrix.ravel()
    if (not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values)
            or not np.isfinite(values).all() or (values < 0).any()
            or (values != np.floor(values)).any() or (values >= 2**53).any()):
        raise ValueError("Nonnegative finite raw integer counts below 2**53 required, not normalized/log values")
    largest = int(values.max()) if len(values) else 0
    if largest * matrix.shape[1] > np.iinfo(np.int64).max:
        raise ValueError("Potential int64 count overflow; no unsafe projection")
    return sparse.csr_matrix(matrix, dtype=np.int64, copy=True)


def _feature_axis(adata, feature_key, strip_versions):
    ids = pd.Index(adata.var_names if feature_key is None else adata.var[feature_key])
    if ids.hasnans or not all(isinstance(x, str) and x for x in ids):
        raise ValueError("Feature IDs must be nonempty strings")
    if strip_versions:
        ids = pd.Index([re.sub(r"^(ENS[A-Z]*G\d+)\.\d+$", r"\1", x) for x in ids])
    if not ids.is_unique:
        raise ValueError("Duplicate input features are not silently merged, including version-stripped collisions")
    return ids


def _resolve_symbols(edges):
    candidates = edges.loc[edges.source_symbol.ne("")].copy()
    if candidates.groupby("source_symbol").source_id.nunique().gt(1).any():
        raise ValueError("Source symbols are not unique stable IDs; use Ensembl features")
    return candidates


def convert_anndata(adata, mapping: MappingTable, *, strategy="one2one", counts_layer="counts",
                    feature_key=None, feature_type="ensembl", output_ids="ensembl",
                    source_ambiguity="error", tie_policy="error", scores=None,
                    min_score=None, strip_ensembl_versions=False, additional_count_layers=()):
    """Map mouse/fish raw counts to human, returning a new AnnData and full decisions.

    one2one: final reciprocal pairs only. sum: sum eligible source paralogs.
    best_sequence/best_expression: one present source feature per human target.
    Sources linked to multiple humans are rejected or explicitly dropped in ALL
    modes; this version never copies one source count to several human genes.
    """
    if strategy not in {"one2one", "sum", "best_sequence", "best_expression", "highest_total_expression"}:
        raise ValueError("Unknown mapping strategy")
    if source_ambiguity not in {"error", "drop"} or tie_policy not in {"error", "drop"}:
        raise ValueError("Ambiguity/tie policies must be error or drop")
    if feature_type not in {"ensembl", "symbol"} or output_ids not in {"ensembl", "symbol"}:
        raise ValueError("Feature/output IDs must be ensembl or symbol")
    if not adata.obs_names.is_unique:
        raise ValueError("Cell barcodes must be unique; no automatic renaming")
    if counts_layer is not None and counts_layer not in adata.layers:
        raise ValueError("Count layer missing; pass counts_layer=None explicitly for raw X")
    matrix = raw_counts(adata.X if counts_layer is None else adata.layers[counts_layer])
    axis = _feature_axis(adata, feature_key, strip_ensembl_versions)
    edges = _canonical(mapping.edges)
    if feature_type == "symbol":
        edges = _resolve_symbols(edges)
    lookup = "source_id" if feature_type == "ensembl" else "source_symbol"
    decisions = edges.loc[edges[lookup].isin(axis)].copy().reset_index(drop=True)
    decisions["input_feature"] = decisions[lookup]
    decisions["input_column"] = axis.get_indexer(decisions[lookup])
    decisions["selected"] = False
    decisions["decision"] = "candidate"
    if strategy == "highest_total_expression":
        largest = int(matrix.data.max()) if matrix.nnz else 0
        if largest * matrix.shape[0] > np.iinfo(np.int64).max:
            raise ValueError("Potential int64 overflow in whole-input expression totals")
        totals = np.asarray(matrix.sum(axis=0)).ravel()
        if (totals >= 2**53).any():
            raise ValueError("Expression totals exceed exact float score representation")
        decisions["total_expression"] = totals[decisions.input_column]
    if strategy == "one2one":
        removed = decisions.human_perspective_type.ne("One2One")
        decisions.loc[removed, "decision"] = "not_final_reciprocal_One2One"
    active = decisions.decision.eq("candidate")
    ambiguous = active & decisions.source_degree.gt(1)
    if ambiguous.any() and source_ambiguity == "error":
        ids = decisions.loc[ambiguous, "source_id"].unique()[:5]
        raise ValueError("Source genes link to multiple human targets; use source_ambiguity='drop' "
                         "or a reviewed revised mapping. Examples: " + ", ".join(ids))
    decisions.loc[ambiguous, "decision"] = "ambiguous_source_dropped_no_count_duplication"
    if scores is not None:
        scores = scores.copy(deep=True)
        if strategy != "best_expression":
            raise ValueError("External scores are supported only for best_expression")
        needed = {"source_id", "target_id", "expression_similarity"}
        if not needed.issubset(scores.columns) or scores.duplicated(["source_id", "target_id"]).any():
            raise ValueError("Expression scores need unique source_id/target_id pairs")
        score_values = pd.to_numeric(scores.expression_similarity, errors="raise")
        if np.isinf(score_values).any() or ((score_values.dropna() < -1) | (score_values.dropna() > 1)).any():
            raise ValueError("Expression scores must be correlations in [-1,1], or NaN")
        scores["expression_similarity"] = score_values
        decisions = decisions.drop(columns=["expression_similarity"], errors="ignore").merge(
            scores[list(needed)], on=["source_id", "target_id"], how="left", validate="one_to_one", sort=False)
    score_key = ("expression_similarity" if strategy == "best_expression" else
                 "total_expression" if strategy == "highest_total_expression" else "sequence_score")
    if strategy == "best_expression" and score_key not in decisions:
        raise ValueError("best_expression requires explicit matched-group expression-similarity scores")
    if min_score is not None and not np.isfinite(min_score):
        raise ValueError("min_score must be finite")
    active = decisions.decision.eq("candidate")
    for target, group in decisions.loc[active].groupby("target_id", sort=True):
        positions = group.index
        if strategy in {"one2one", "sum"} or len(group) == 1:
            winners = positions
        else:
            if score_key not in group or group[score_key].isna().any():
                if tie_policy == "error":
                    raise ValueError("Missing comparison scores for human target " + target)
                decisions.loc[positions, "decision"] = "missing_comparison_score_target_dropped"
                continue
            values = group[score_key].to_numpy(dtype=float)
            winners = positions[np.isclose(values, values.max(), rtol=0, atol=1e-12)]
            if len(winners) != 1:
                if tie_policy == "error":
                    raise ValueError("Tied best scores for human target " + target)
                decisions.loc[positions, "decision"] = "tied_score_target_dropped"
                continue
            decisions.loc[positions.difference(winners), "decision"] = "lower_representative_score"
        if min_score is not None:
            if strategy in {"one2one", "sum"}:
                raise ValueError("min_score applies only to representative-selection strategies")
            values = decisions.loc[winners, score_key] if score_key in decisions else pd.Series(np.nan, index=winners)
            if values.isna().any() or (values < min_score).any():
                decisions.loc[positions, "decision"] = "below_min_score_target_dropped"
                continue
        decisions.loc[winners, "selected"] = True
        decisions.loc[winners, "decision"] = "selected"
    selected = decisions.loc[decisions.selected].copy()
    if selected.empty:
        raise ValueError("No genes survive the explicit mapping/ambiguity policies")
    if selected.input_column.duplicated().any():
        raise ValueError("A source feature would be counted more than once")
    targets = pd.Index(sorted(selected.target_id.unique()))
    projection = sparse.csr_matrix((np.ones(len(selected), dtype=np.int64),
        (selected.input_column, targets.get_indexer(selected.target_id))), shape=(len(axis), len(targets)))
    mapped = matrix @ projection
    selected_input_counts = np.asarray(matrix[:, selected.input_column].sum(axis=1)).ravel()
    if not np.array_equal(np.asarray(mapped.sum(axis=1)).ravel(), selected_input_counts):
        raise ValueError("Count accounting failed; mapped mass must equal selected input mass")
    metadata = selected.groupby("target_id", sort=True).agg(target_symbol=("target_symbol", "first"),
        source_ids=("source_id", lambda s: ";".join(sorted(s))),
        source_features=("input_feature", lambda s: ";".join(sorted(s))), n_source_features=("source_id", "size"))
    metadata = metadata.loc[targets]
    metadata.index.name = "human_ensembl_id"
    metadata["human_ensembl_id"] = metadata.index
    if output_ids == "symbol":
        symbols = pd.Index(metadata.target_symbol)
        if (symbols == "").any() or not symbols.is_unique:
            raise ValueError("Human symbols are empty/nonunique; keep stable-ID output")
        metadata.index = symbols
    result = ad.AnnData(X=mapped, obs=adata.obs.copy(deep=True), var=metadata,
        obsm={k: deepcopy(v) for k, v in adata.obsm.items()},
        obsp={k: deepcopy(v) for k, v in adata.obsp.items()})
    result.layers["counts"] = mapped.copy()
    for layer in additional_count_layers:
        if layer == "counts" or layer == counts_layer:
            continue
        if layer not in adata.layers:
            raise ValueError("Additional raw count layer is missing: " + layer)
        result.layers[layer] = raw_counts(adata.layers[layer]) @ projection
    library = np.asarray(matrix.sum(axis=1)).ravel()
    reserved = "heart_orthomap_input_library_size"
    if reserved in adata.obs and not np.array_equal(adata.obs[reserved].to_numpy(), library):
        raise ValueError("Reserved input-library column already exists with different values")
    result.obs[reserved] = library
    mapped_sources = set(selected.input_feature)
    report = {"strategy": strategy, "source_species": mapping.source_species, "target_species": "human",
        "counts_layer": "X" if counts_layer is None else counts_layer,
        "feature_type": feature_type, "feature_key": feature_key or "var_names", "output_ids": output_ids,
        "ensembl_version_stripping": strip_ensembl_versions,
        "min_score": "not_requested" if min_score is None else float(min_score),
        "source_ambiguity": source_ambiguity, "tie_policy": tie_policy,
        "input_cells": adata.n_obs, "input_features": adata.n_vars,
        "selected_source_features": len(mapped_sources), "output_human_features": len(targets),
        "aggregated_human_features": int(metadata.n_source_features.gt(1).sum()),
        "unmapped_or_excluded_input_features": len(axis) - len(mapped_sources),
        "cell_barcodes_and_obs_preserved": True, "source_counts_modified": False,
        "source_count_duplication": False, "selected_count_mass_preserved_per_cell": True,
        "mapping": deepcopy(mapping.provenance),
        "mapping_rows_quarantined_for_identifier_review": len(mapping.issues),
        "dropped_feature_axis_objects": [k for k, yes in (("raw", adata.raw is not None),
            ("varm", bool(adata.varm)), ("varp", bool(adata.varp))) if yes],
        "dropped_layers": sorted(set(adata.layers) - set(additional_count_layers) - {counts_layer, "counts"}),
        "original_uns_not_copied": sorted(adata.uns.keys()),
        "normalization": "No normalization performed. Use original full-native library sizes for comparisons.",
        "obsm_obsp": "Preserved on identical cell axes; native embeddings/graphs are not recomputed",
        "unmapped_input_features": [x for x in axis if x not in mapped_sources]}
    if scores is not None:
        normalized_scores = scores[["source_id", "target_id", "expression_similarity"]].sort_values(
            ["source_id", "target_id"]).to_csv(index=False, sep="\t", lineterminator="\n")
        report["expression_scores_sha256"] = hashlib.sha256(normalized_scores.encode()).hexdigest()
        report["expression_score_provenance"] = deepcopy(scores.attrs.get("heart_orthomap_scoring", {
            "method": "Externally supplied pairwise correlations; caller must retain training/reference provenance"}))
    result.uns["heart_orthomap"] = deepcopy(report)
    return ConversionResult(result, decisions, report)
