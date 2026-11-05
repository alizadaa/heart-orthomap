"""Load frozen fixtures and create fresh inputs for every independent check."""

import copy
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from validation_contract import Request

ROOT = Path(__file__).resolve().parent


def load_json(name):
    with (ROOT / "fixtures" / name).open() as handle:
        return json.load(handle)


def core_request(method="sum", fmt="csr"):
    fixture = load_json("core.json")
    matrix = sparse.coo_matrix(np.asarray(fixture["counts"], dtype=np.int64)).asformat(fmt)
    matrix.sum_duplicates()
    symbols = dict(zip(fixture["source_ids"], fixture["source_symbols"]))
    for edge in fixture["edges"]:
        edge["source_symbol"] = symbols.get(edge["source_id"], edge["source_id"])
    request = Request(matrix, tuple(fixture["source_ids"]), tuple(fixture["edges"]), method)
    if method == "best_expression":
        request.scores = [{"source_id": edge["source_id"], "target_id": edge["human_id"],
                           "expression_similarity": 0.9 if edge["source_id"] == "S021" else 0.1}
                          for edge in fixture["edges"]]
        unique = {(row["source_id"], row["target_id"]): row for row in request.scores}
        request.scores = list(unique.values())
    return request


def pair_edge(human, source, score=None, kind=None):
    edge = {"human_id": human, "source_id": source, "sequence_score": score}
    if kind is not None:
        edge["orthology_type"] = kind
    return edge


def small_request(values, source_ids=None, edges=None, method="sum", dtype=np.int64):
    matrix = sparse.csr_matrix(np.asarray(values, dtype=dtype))
    sources = tuple(source_ids if source_ids is not None else (
        "S%03d" % index for index in range(matrix.shape[1])
    ))
    if edges is None:
        edges = [pair_edge("H", source, 90) for source in sources]
    edges = copy.deepcopy(edges)
    source_humans, human_sources = {}, {}
    for edge in edges:
        source_humans.setdefault(edge["source_id"], set()).add(edge["human_id"])
        human_sources.setdefault(edge["human_id"], set()).add(edge["source_id"])
    for edge in edges:
        source_degree = len(source_humans[edge["source_id"]])
        target_degree = len(human_sources[edge["human_id"]])
        kind = ("one2one" if source_degree == target_degree == 1 else
                "one2many" if source_degree == 1 else "many2one" if target_degree == 1 else "many2many")
        edge.setdefault("orthology_type", kind)
    return Request(matrix, sources, tuple(edges), method)


def selection_requests():
    for case in load_json("selection.json")["cases"]:
        edges = [pair_edge("H", source, float("nan") if score == "NaN" else score)
                 for source, score in zip(case["source_ids"], case["scores"])]
        for method, key in (("best_sequence", "sequence"), ("best_expression", "expression")):
            yield case, key, small_request(case["counts"], case["source_ids"], edges, method)


def metadata_request(method="best_expression"):
    import anndata as ad
    import pandas as pd

    fixture = load_json("core.json")
    request = core_request(method)
    transformed = request.X.astype(np.float64).toarray()
    transformed[:, 1] = [0.01, 0.01, 0.01]
    transformed[:, 2] = [100.0, 100.0, 100.0]
    obs = pd.DataFrame({
        "stage": pd.Categorical(["late", "early", "late"], categories=["early", "late"], ordered=True),
        "donor": ["donor_2", "donor_1", "donor_2"],
        "quality": [0.9, 0.4, np.nan],
    }, index=fixture["obs_ids"])
    var = pd.DataFrame({
        "source_id": fixture["source_ids"],
        "gene_symbol": fixture["source_symbols"],
        "native_annotation": list(range(len(fixture["source_ids"]))),
    }, index=fixture["source_ids"])
    data = ad.AnnData(sparse.csr_matrix(transformed), obs=obs, var=var)
    data.layers["counts"] = request.X.copy()
    data.layers["alt_counts"] = sparse.csr_matrix(np.asarray([
        [5, 8, 3, 1, 2, 3, 4, 0],
        [0, 1, 100, 2, 3, 4, 5, 0],
        [7, 2, 4, 3, 4, 5, 6, 0],
    ], dtype=np.int64))
    data.obsm["X_embedding"] = np.asarray([[91, 2], [-4, 3], [0, 5]], dtype=float)
    data.obsp["connectivities"] = sparse.csr_matrix([[0, 2, 0], [2, 0, 3], [0, 3, 0]])
    data.uns["provenance"] = {"native_species": "synthetic", "nested": {"batches": [2, 1]}}
    data.varm["native_loadings"] = np.arange(data.n_vars * 2).reshape(data.n_vars, 2)
    data.varp["native_gene_graph"] = sparse.eye(data.n_vars, format="csr")
    data.raw = data.copy()
    request.adata = data
    request.count_layer = "counts"
    request.additional_count_layers = ("alt_counts",)
    return request


def snapshot(request):
    """Capture meaningful input state, excluding harmless sparse cache flags."""
    def matrix_state(matrix):
        state = [matrix.format, matrix.shape, matrix.dtype.str]
        for name in ("data", "indices", "indptr", "row", "col"):
            if hasattr(matrix, name):
                array = getattr(matrix, name)
                state.append((name, array.dtype.str, array.shape, array.tobytes()))
        return state

    state = {"X": matrix_state(request.X), "ids": copy.deepcopy(request.source_ids),
             "edges": copy.deepcopy(request.edges), "scores": copy.deepcopy(request.scores)}
    if request.adata is not None:
        data = request.adata
        state["adata"] = {
            "X": matrix_state(data.X), "obs": data.obs.copy(deep=True),
            "var": data.var.copy(deep=True), "uns": copy.deepcopy(data.uns),
            "layers": {key: matrix_state(data.layers[key]) for key in data.layers},
            "obsm": {key: data.obsm[key].copy() for key in data.obsm},
            "obsp": {key: matrix_state(data.obsp[key]) for key in data.obsp},
            "varm": {key: data.varm[key].copy() for key in data.varm},
            "varp": {key: matrix_state(data.varp[key]) for key in data.varp},
            "raw_X": matrix_state(data.raw.X), "raw_var": data.raw.var.copy(deep=True),
        }
    return state
