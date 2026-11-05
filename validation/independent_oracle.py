"""Small exact-integer oracle, independent of sparse mapping products."""

import copy
import math
import numbers

import numpy as np
from scipy import sparse

from validation_contract import Result

INT64_MAX = 2**63 - 1
FLOAT_SAFE_MAX = 2**53
METHODS = ("one2one", "sum", "best_sequence", "best_expression", "highest_total_expression")


def valid_id(value):
    return isinstance(value, str) and bool(value.strip())


def integer_entries(matrix, source_ids):
    if not sparse.issparse(matrix) or matrix.ndim != 2:
        raise ValueError("Expected a two-dimensional sparse count matrix")
    if matrix.shape[1] != len(source_ids):
        raise ValueError("Feature IDs do not match matrix columns")
    if not all(valid_id(value) for value in source_ids):
        raise ValueError("Missing or invalid stable source ID")
    if len(set(source_ids)) != len(source_ids):
        raise ValueError("Duplicate stable source IDs")
    if matrix.dtype.kind not in "iuf":
        raise ValueError("Raw counts must be real integer-valued numbers")
    if not getattr(matrix, "has_canonical_format", True):
        raise ValueError("Noncanonical sparse counts are explicitly rejected")
    coo = matrix.tocoo(copy=True)
    entries = {}
    for row, col, raw in zip(coo.row, coo.col, coo.data):
        if isinstance(raw, numbers.Integral):
            value = int(raw)
        else:
            value_float = float(raw)
            if not math.isfinite(value_float) or not value_float.is_integer():
                raise ValueError("Nonfinite or fractional raw count")
            if value_float >= FLOAT_SAFE_MAX:
                raise ValueError("Floating count exceeds safe exact-integer range")
            value = int(value_float)
        if value < 0:
            raise ValueError("Negative raw count")
        if value >= FLOAT_SAFE_MAX:
            raise ValueError("Count entries must be below 2**53")
        if value > INT64_MAX:
            raise OverflowError("Raw count exceeds signed int64 range")
        key = (int(row), int(col))
        total = entries.get(key, 0) + value
        if total > INT64_MAX:
            raise OverflowError("Duplicate-coordinate sum exceeds signed int64 range")
        if total:
            entries[key] = total
    return entries


def score_value(raw):
    if raw is None:
        return None
    value = float(raw)
    if math.isnan(value):
        return None
    if not math.isfinite(value):
        raise ValueError("Infinite sequence score")
    return value


def selected_plan(request, entries):
    if request.method not in METHODS:
        raise ValueError("Unknown aggregation method")
    if request.ambiguity not in ("drop", "error"):
        raise ValueError("Unknown ambiguity policy")
    present = set(request.source_ids)
    pair_scores = {}
    pair_types = {}
    for edge in request.edges:
        human, source = edge.get("human_id"), edge.get("source_id")
        if not valid_id(human) or not valid_id(source):
            raise ValueError("Missing or invalid mapping endpoint ID")
        score = score_value(edge.get("sequence_score"))
        if source not in present:
            continue
        pair = (human, source)
        previous = pair_scores.get(pair)
        if pair in pair_scores and score is not None and previous is not None and score != previous:
            raise ValueError("Conflicting repeated-pair sequence scores")
        if pair not in pair_scores or score is not None:
            pair_scores[pair] = score
        pair_types[pair] = edge.get("orthology_type", "one2many")
    source_humans = {}
    for human, source in pair_scores:
        source_humans.setdefault(source, set()).add(human)
    ambiguous = sorted(source for source, humans in source_humans.items() if len(humans) > 1)
    if request.method == "one2one":
        eligible_pairs = {pair for pair in pair_scores if pair_types[pair].lower() == "one2one"}
        ambiguous = [source for source in ambiguous if any(pair[1] == source for pair in eligible_pairs)]
    else:
        eligible_pairs = set(pair_scores)
    if ambiguous and request.ambiguity == "error":
        raise ValueError("Ambiguous source mappings: " + ", ".join(ambiguous))
    totals = {source: 0 for source in request.source_ids}
    for (_, col), count in entries.items():
        totals[request.source_ids[col]] += count
    candidates = {}
    for human, source in pair_scores:
        if (human, source) in eligible_pairs and len(source_humans[source]) == 1:
            candidates.setdefault(human, []).append(source)
    expression = {}
    if request.scores is not None:
        for row in request.scores:
            pair = (row["target_id"], row["source_id"])
            if pair in expression:
                raise ValueError("Duplicate expression-score pairs")
            expression[pair] = score_value(row.get("expression_similarity"))
    if request.method == "best_expression" and request.scores is None:
        raise ValueError("Explicit expression similarity is required")
    selected = {}
    for human in sorted(candidates):
        sources = sorted(candidates[human])
        if request.method in ("best_sequence", "best_expression", "highest_total_expression"):
            values = {source: (pair_scores[(human, source)] if request.method == "best_sequence"
                               else totals[source] if request.method == "highest_total_expression"
                               else expression.get((human, source))) for source in sources}
            if request.method == "highest_total_expression" and any(value >= 2**53 for value in totals.values()):
                raise ValueError("Expression totals exceed safe score range")
            if len(sources) > 1:
                missing = any(value is None for value in values.values())
                winners = [] if missing else [source for source in sources
                    if abs(values[source] - max(values.values())) <= 1e-12]
                if missing or len(winners) != 1:
                    if request.tie_policy == "error":
                        raise ValueError("Missing or tied comparison scores")
                    continue
                sources = winners
            if request.min_score is not None:
                if any(values[source] is None or values[source] < request.min_score for source in sources):
                    continue
        selected[human] = tuple(sources)
    if not selected:
        raise ValueError("No genes survive mapping policies")
    return selected, {
        "ambiguous_source_ids": ambiguous,
        "unmatched_source_ids": sorted(present - set(source_humans)),
    }


def apply_plan(matrix, source_ids, selected):
    entries = integer_entries(matrix, source_ids)
    humans = tuple(selected)
    source_destination = {
        source: column
        for column, human in enumerate(humans)
        for source in selected[human]
    }
    output = {}
    for (row, col), value in entries.items():
        source = source_ids[col]
        if source not in source_destination:
            continue
        key = (row, source_destination[source])
        total = output.get(key, 0) + value
        if total > INT64_MAX:
            raise OverflowError("Mapped sum exceeds signed int64 range")
        output[key] = total
    coordinates = sorted(output)
    return sparse.coo_matrix((
        np.asarray([output[key] for key in coordinates], dtype=np.int64),
        ([key[0] for key in coordinates], [key[1] for key in coordinates]),
    ), shape=(matrix.shape[0], len(humans)), dtype=np.int64).tocsr()


def map_counts(request):
    entries = integer_entries(request.X, request.source_ids)
    selected, provenance = selected_plan(request, entries)
    mapped = apply_plan(request.X, request.source_ids, selected)
    result = Result(mapped, tuple(selected), selected_sources=selected, provenance=provenance)
    if request.adata is not None:
        import anndata as ad
        import pandas as pd

        original = request.adata
        if original.shape != request.X.shape:
            raise ValueError("AnnData and mapping matrix shapes disagree")
        if request.count_layer not in original.layers:
            raise ValueError("Requested count layer is absent")
        output = ad.AnnData(
            X=mapped.copy(),
            obs=original.obs.copy(deep=True),
            var=pd.DataFrame({"human_id": list(selected)}, index=list(selected)),
        )
        for key in original.obsm:
            output.obsm[key] = original.obsm[key].copy()
        for key in original.obsp:
            output.obsp[key] = original.obsp[key].copy()
        for key in (request.count_layer,) + request.additional_count_layers:
            if key not in original.layers:
                raise ValueError("Requested additional count layer is absent")
            output.layers[key] = apply_plan(original.layers[key], request.source_ids, selected)
        result.adata = output
    return result
