"""Explicit matched-group expression profiles, never correlations of unpaired cells."""
import numpy as np
import pandas as pd
from scipy import sparse

from .convert import raw_counts, _feature_axis


def _profiles(adata, groups, layer, feature_key, min_cells):
    matrix = raw_counts(adata.X if layer is None else adata.layers[layer])
    ids = _feature_axis(adata, feature_key, False)
    if isinstance(groups, pd.Series):
        if not groups.index.equals(adata.obs_names):
            raise ValueError("Group Series must have the exact ordered cell index")
        values = groups.to_numpy()
    else:
        values = np.asarray(groups)
    if values.shape != (adata.n_obs,) or pd.isna(values).any() or not all(isinstance(x, str) and x for x in values):
        raise ValueError("Supply one nonempty matched-group label per cell")
    library = np.asarray(matrix.sum(axis=1)).ravel().astype(float)
    if (library <= 0).any():
        raise ValueError("Expression scoring requires positive full-native input libraries")
    labels, counts = np.unique(values, return_counts=True)
    labels = labels[counts >= min_cells]
    count_lookup = {label: int(n) for label, n in zip(*np.unique(values, return_counts=True))}
    supported = np.isin(values, labels)
    positions = pd.Index(labels).get_indexer(values[supported])
    weights = sparse.csr_matrix((1 / np.array([count_lookup[label] for label in values[supported]]),
                                 (positions, np.flatnonzero(supported))), shape=(len(labels), adata.n_obs))
    normalized = sparse.diags(10000 / library) @ matrix
    profiles = (weights @ normalized).toarray()
    return pd.DataFrame(np.log1p(profiles), index=labels, columns=ids), {label: count_lookup[label] for label in labels}


def expression_similarity_scores(source, reference, mapping, *, source_groups, reference_groups,
                                 source_counts_layer="counts", reference_counts_layer="counts",
                                 source_feature_key=None, reference_feature_key=None,
                                 min_cells=5, min_groups=3):
    """Pearson r over log1p(mean native CP10K) in explicitly matched groups.

    Both objects must retain the full native gene axis used for library sizes.
    Group names must already represent biological correspondences, e.g. matched
    cell-type/stage strata, not arbitrary species-specific cluster numbers.
    This helper selects no genes and must not be trained on benchmark holdouts.
    """
    if isinstance(min_cells, bool) or not isinstance(min_cells, int) or min_cells < 1:
        raise ValueError("min_cells must be a positive integer")
    if isinstance(min_groups, bool) or not isinstance(min_groups, int) or min_groups < 3:
        raise ValueError("At least three matched groups are required")
    left, left_n = _profiles(source, source_groups, source_counts_layer, source_feature_key, min_cells)
    right, right_n = _profiles(reference, reference_groups, reference_counts_layer, reference_feature_key, min_cells)
    common = left.index.intersection(right.index).sort_values()
    if len(common) < min_groups:
        raise ValueError("Insufficient adequately supported matched expression groups")
    rows = []
    for edge in mapping.edges.itertuples(index=False):
        if edge.source_id not in left or edge.target_id not in right:
            continue
        x, y = left.loc[common, edge.source_id].to_numpy(), right.loc[common, edge.target_id].to_numpy()
        constant = np.ptp(x) == 0 or np.ptp(y) == 0
        x, y = x - x.mean(), y - y.mean()
        denominator = np.linalg.norm(x) * np.linalg.norm(y)
        score = float(np.clip(x @ y / denominator, -1, 1)) if denominator > 0 and not constant else np.nan
        rows.append({"source_id": edge.source_id, "target_id": edge.target_id,
                     "expression_similarity": score, "matched_groups": len(common)})
    result = pd.DataFrame(rows, columns=["source_id", "target_id", "expression_similarity", "matched_groups"])
    result.attrs["heart_orthomap_scoring"] = {
        "method": "Pearson r of log1p(mean native CP10K) across matched groups",
        "matched_groups": common.tolist(), "source_cells_per_group": {g: left_n[g] for g in common},
        "reference_cells_per_group": {g: right_n[g] for g in common},
        "group_weighting": "Equal groups; pooled cells within each group; no donor balancing",
        "library_scope": "All features in the supplied native objects; supply full native axes",
        "benchmark_caveat": "Use independent training/reference cells, not held-out evaluation labels"}
    return result
