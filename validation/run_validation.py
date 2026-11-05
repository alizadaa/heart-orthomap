#!/usr/bin/env python3
"""Independent executable checks. Writes only beneath this script's directory."""

import argparse
import copy
from dataclasses import replace
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import time
import traceback
import warnings

sys.dont_write_bytecode = True
if not __debug__:
    raise RuntimeError("Assertions are required; run without -O or PYTHONOPTIMIZE")

import numpy as np
import pandas as pd
from scipy import sparse

from fixture_support import ROOT, core_request, load_json, metadata_request, pair_edge, small_request, snapshot
from validation_contract import AdapterNotReady, Request, Result
import independent_oracle


def equal(left, right):
    if isinstance(left, pd.DataFrame):
        pd.testing.assert_frame_equal(left, right)
    elif isinstance(left, np.ndarray):
        np.testing.assert_array_equal(left, right)
    elif isinstance(left, dict):
        assert left.keys() == right.keys(), (left.keys(), right.keys())
        for key in left:
            equal(left[key], right[key])
    elif isinstance(left, (list, tuple)):
        assert len(left) == len(right)
        for a, b in zip(left, right):
            equal(a, b)
    elif isinstance(left, float) and np.isnan(left):
        assert isinstance(right, float) and np.isnan(right)
    else:
        assert left == right, (left, right)


def matrix_values(matrix):
    assert sparse.issparse(matrix), "Output count matrix must remain sparse"
    assert matrix.dtype.kind in "iu", "Output raw counts must have integer dtype"
    values = {}
    coo = matrix.tocoo(copy=True)
    for row, col, value in zip(coo.row, coo.col, coo.data):
        assert int(value) >= 0, "Negative/wrapped output count"
        key = (int(row), int(col))
        values[key] = values.get(key, 0) + int(value)
    return {key: value for key, value in values.items() if value}


def assert_matrix(matrix, expected):
    if not sparse.issparse(expected):
        expected = sparse.csr_matrix(np.asarray(expected, dtype=np.int64))
    assert matrix.shape == expected.shape, (matrix.shape, expected.shape)
    assert matrix_values(matrix) == matrix_values(expected), (matrix_values(matrix), matrix_values(expected))


def assert_result(result, humans, expected):
    assert len(result.human_ids) == len(set(result.human_ids)), "Duplicate output human IDs"
    assert set(result.human_ids) == set(humans), (result.human_ids, humans)
    order = [result.human_ids.index(human) for human in humans]
    assert_matrix(result.X.tocsr()[:, order], expected)


def reject(function, request):
    before = snapshot(request)
    try:
        function(request)
    except (ValueError, OverflowError) as error:
        assert str(error), "Rejection needs an informative error"
    else:
        raise AssertionError("Invalid/unsupported input silently accepted")
    finally:
        equal(before, snapshot(request))


class Suite:
    def __init__(self, mapper, package=None):
        self.mapper = mapper
        self.package = package
        self.results = []

    def case(self, name, check, priority="P0"):
        start = time.monotonic()
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                check()
            item = {"id": name, "status": "passed", "priority": priority,
                    "warnings": [str(w.message) for w in caught]}
        except Exception as error:
            item = {"id": name, "status": "failed", "priority": priority,
                    "error_type": type(error).__name__, "error": str(error),
                    "traceback": traceback.format_exc()}
        item["seconds"] = round(time.monotonic() - start, 4)
        self.results.append(item)
        print(item["status"].upper(), name, item.get("error", ""), flush=True)

    def checked(self, request):
        before = snapshot(request)
        try:
            result = self.mapper(request)
        finally:
            equal(before, snapshot(request))
        assert result.X.shape[0] == request.X.shape[0], "Cell axis changed"
        input_totals = {}
        for (row, _), value in matrix_values(request.X.astype(np.int64)).items():
            input_totals[row] = input_totals.get(row, 0) + value
        output_totals = {}
        for (row, _), value in matrix_values(result.X).items():
            output_totals[row] = output_totals.get(row, 0) + value
        assert all(value <= input_totals.get(row, 0) for row, value in output_totals.items()), "Count inflation"
        if result.selected_sources is not None:
            sources = [source for group in result.selected_sources.values() for source in group]
            assert len(sources) == len(set(sources)), "Source selected for multiple humans"
        return result

    def run(self):
        golden = load_json("golden.json")["core"]
        for method in ("sum", "best_sequence", "best_expression", "highest_total_expression"):
            for fmt in ("csr", "csc", "coo"):
                def check(method=method, fmt=fmt):
                    request = core_request(method, fmt)
                    if self.package is not None and fmt == "coo":
                        import anndata as ad
                        try:
                            ad.AnnData(X=request.X)
                        except ValueError as error:
                            if "Only CSR and CSC" not in str(error):
                                raise
                            # New AnnData rejects COO before conversion is called.
                            reject(self.mapper, request)
                            return
                    result = self.checked(request)
                    assert_result(result, golden["human_ids"], golden[method])
                    equal(result.selected_sources, {k: tuple(v) for k, v in golden["selected_" + method].items()})
                self.case("T01_T04_core_%s_%s" % (method, fmt), check)
            for seed in range(5):
                def check(method=method, seed=seed):
                    request = core_request(method)
                    rng = np.random.default_rng(seed)
                    source_order = rng.permutation(len(request.source_ids))
                    obs_order = rng.permutation(request.X.shape[0])
                    request.X = request.X[obs_order, :][:, source_order]
                    request.X.sort_indices()
                    request.source_ids = tuple(request.source_ids[index] for index in source_order)
                    edges = list(request.edges) * (seed + 2)
                    rng.shuffle(edges)
                    request.edges = tuple(edges)
                    result = self.checked(request)
                    assert_result(result, golden["human_ids"], np.asarray(golden[method])[obs_order, :])
                self.case("T02_T14_permuted_repeated_%s_seed%d" % (method, seed), check)
        self.case("T03_ambiguous_error", lambda: reject(self.mapper, replace(core_request(), ambiguity="error")))

        def strict():
            result = self.checked(core_request("one2one"))
            assert_result(result, ["H010", "H060"], [[2, 0], [0, 0], [4, 0]])
        self.case("T01_strict_one2one_direction", strict)

        for fixture in load_json("selection.json")["cases"]:
            for policy in ("error", "drop"):
                def check(fixture=fixture, policy=policy):
                    edges = [pair_edge("H", source, None if score == "NaN" else score)
                             for source, score in zip(fixture["source_ids"], fixture["scores"])]
                    # A separate sentinel survives a target-level drop, making zero survivors distinct.
                    sources = fixture["source_ids"] + ["S_SENTINEL"]
                    counts = [row + [3] for row in fixture["counts"]]
                    edges.append(pair_edge("H_SENTINEL", "S_SENTINEL", 90, "one2one"))
                    request = small_request(counts, sources, edges, "best_sequence")
                    request.tie_policy = policy
                    disputed = fixture["name"] != "global_source_not_per_cell_max"
                    if fixture["name"] == "expression_tie":
                        request.method = "best_expression"
                        request.scores = [{"source_id": source, "target_id": "H", "expression_similarity": 0.5}
                                          for source in fixture["source_ids"]]
                        request.scores.append({"source_id": "S_SENTINEL", "target_id": "H_SENTINEL", "expression_similarity": 0.8})
                    if disputed and policy == "error":
                        reject(self.mapper, request)
                    elif disputed:
                        assert_result(self.checked(request), ["H_SENTINEL"], [[3], [3]])
                    else:
                        result = self.checked(request)
                        assert_result(result, ["H", "H_SENTINEL"], [[1, 3], [7, 3]])
                self.case("T05_T06_%s_%s" % (fixture["name"], policy), check, "P1")

        def explicit_expression():
            request = small_request([[100, 1], [50, 2]], ["S_A", "S_B"], method="best_expression")
            request.scores = [{"source_id": "S_A", "target_id": "H", "expression_similarity": -0.2},
                              {"source_id": "S_B", "target_id": "H", "expression_similarity": 0.8}]
            assert_result(self.checked(request), ["H"], [[1], [2]])
        self.case("T07_expression_similarity_not_abundance", explicit_expression)
        self.case("T07_expression_scores_required", lambda: reject(self.mapper, small_request([[1, 2]], method="best_expression")), "P1")

        for label, value, dtype, accepted in (
            ("int32_promotion", 2**30, np.int32, True),
            ("safe_large_int", 2**52 + 1, np.int64, True),
            ("whole_float", 2.0, np.float64, True),
            ("negative", -1, np.int64, False),
            ("fractional", 0.5, np.float64, False),
            ("nan", np.nan, np.float64, False),
            ("positive_inf", np.inf, np.float64, False),
            ("negative_inf", -np.inf, np.float64, False),
            ("int_above_2p53", 2**53 + 1, np.int64, False),
            ("float_at_2p53", float(2**53), np.float64, False),
            ("uint64_boundary", 2**63, np.uint64, False),
            ("int64_overflow", 2**62, np.int64, False),
            ("boolean", True, np.bool_, False),
            ("complex", 1 + 2j, np.complex128, False),
        ):
            def check(value=value, dtype=dtype, accepted=accepted):
                request = small_request([[value, value]], dtype=dtype)
                if accepted:
                    assert_result(self.checked(request), ["H"], [[2 * int(value)]])
                else:
                    reject(self.mapper, request)
            self.case("T09_T12_numeric_" + label, check)

        for value in (np.nan, np.inf, -1, 0.25):
            def check(value=value):
                request = small_request([[2, value]], ["S_GOOD", "S_UNMATCHED"],
                                        [pair_edge("H", "S_GOOD", 90)], dtype=np.float64)
                reject(self.mapper, request)
            self.case("T11_invalid_unmatched_" + str(value), check)

        def duplicate_sparse():
            request = small_request([[1]])
            request.X = sparse.csr_matrix((np.asarray([3, -2]), np.asarray([0, 0]), np.asarray([0, 2])), shape=(1, 1))
            reject(self.mapper, request)
        self.case("T13_noncanonical_cancellation_rejected", duplicate_sparse)

        for kind in ("no_edges", "no_matching_sources", "ambiguous_only", "zero_sources", "zero_cells"):
            def check(kind=kind):
                if kind == "no_edges":
                    request = small_request([[1]], edges=[])
                elif kind == "no_matching_sources":
                    request = small_request([[1]], edges=[pair_edge("H", "S_ABSENT")])
                elif kind == "ambiguous_only":
                    request = small_request([[1]], ["S"], [pair_edge("H1", "S"), pair_edge("H2", "S")])
                elif kind == "zero_sources":
                    request = Request(sparse.csr_matrix((2, 0), dtype=np.int64), (), (pair_edge("H", "S"),))
                else:
                    request = Request(sparse.csr_matrix((0, 1), dtype=np.int64), ("S",), (pair_edge("H", "S", kind="one2one"),))
                    assert_result(self.checked(request), ["H"], sparse.csr_matrix((0, 1), dtype=np.int64))
                    return
                reject(self.mapper, request)
            self.case("T17_T18_" + kind, check, "P1")

        for label, ids in (("duplicate", ("S", "S")), ("missing", ("S", None)), ("empty", ("S", ""))):
            self.case("T19_source_id_" + label, lambda ids=ids: reject(self.mapper, small_request([[1, 2]], ids)), "P1")
        self.case("T19_shape_mismatch", lambda: reject(self.mapper, small_request([[1, 2]], ["S"])), "P1")

        def conflicting_pair():
            request = small_request([[1]], ["S"], [pair_edge("H", "S", 70), pair_edge("H", "S", 80)])
            reject(self.mapper, request)
        self.case("T20_conflicting_pair_scores", conflicting_pair, "P1")

        for method in ("sum", "best_sequence", "best_expression", "highest_total_expression"):
            def metadata(method=method):
                request = metadata_request(method)
                result = self.checked(request)
                assert_result(result, golden["human_ids"], golden[method])
                output = result.adata
                assert output is not None
                pd.testing.assert_frame_equal(output.obs[request.adata.obs.columns], request.adata.obs)
                equal(output.obsm["X_embedding"], request.adata.obsm["X_embedding"])
                assert_matrix(output.obsp["connectivities"], request.adata.obsp["connectivities"])
                assert_matrix(output.layers["counts"], output.X)
                assert_matrix(output.layers["alt_counts"], load_json("golden.json")["metadata_alt"][method])
                assert output.raw is None and not output.varm and not output.varp, "Stale source feature-axis state"
                assert "provenance" not in output.uns, "Native uns must be dropped/reported under this API"
                if self.package:
                    report = result.provenance
                    assert set(report["dropped_feature_axis_objects"]) == {"raw", "varm", "varp"}
                    assert report["original_uns_not_copied"] == ["provenance"]
                    assert report["source_count_duplication"] is False
                    np.testing.assert_array_equal(output.obs.heart_orthomap_input_library_size, golden["sum_input_totals"])
                    assert report["selected_source_features"] == (4 if method == "sum" else 3)
                before = snapshot(request)
                output.X.data[0] += 1
                equal(before, snapshot(request))
            self.case("T15_T16_metadata_layers_" + method, metadata)

        if self.package:
            self.package_checks()
            from extended_checks import run_extended
            run_extended(self, assert_result, assert_matrix, reject, equal)
        return self.results

    def package_checks(self):
        package = self.package

        def symbol_output():
            request = small_request([[4, 7]], ["S1", "S2"],
                                    [pair_edge("H1", "S1", 90, "one2one"), pair_edge("H2", "S2", 90, "one2one")])
            request.edges = tuple(dict(edge, human_symbol="SYMBOL_" + edge["human_id"]) for edge in request.edges)
            request.output_ids = "symbol"
            result = self.checked(request)
            assert_result(result, ["H1", "H2"], [[4, 7]])
            assert tuple(result.adata.var_names) == ("SYMBOL_H1", "SYMBOL_H2")
        self.case("T08_unique_symbol_output", symbol_output)

        def duplicate_output_symbols():
            request = core_request()
            request.output_ids = "symbol"
            reject(self.mapper, request)
        self.case("T08_duplicate_output_symbols_rejected", duplicate_output_symbols)

        def version_collision():
            request = small_request([[4, 7]], ["ENSMUSG00000000001.1", "ENSMUSG00000000001.2"],
                                    [pair_edge("H", "ENSMUSG00000000001", 90)])
            request.strip_ensembl_versions = True
            reject(self.mapper, request)
        self.case("T19_version_collision", version_collision, "P1")

        def feature_column():
            request = metadata_request("sum")
            request.adata.var_names = ["display_%d" % i for i in range(request.X.shape[1])]
            request.feature_key = "source_id"
            result = self.checked(request)
            golden = load_json("golden.json")["core"]
            assert_result(result, golden["human_ids"], golden["sum"])
        self.case("T14_feature_key_not_var_names", feature_column)

        def roundtrip(method):
            import anndata as ad

            request = metadata_request(method)
            result = self.checked(request)
            golden = load_json("golden.json")["core"]
            if method == "one2one":
                assert_result(result, ["H010", "H060"], [[2, 0], [0, 0], [4, 0]])
            else:
                assert_result(result, golden["human_ids"], golden[method])
            path = ROOT / "reports" / ("roundtrip_" + method + ".h5ad")
            result.adata.write_h5ad(path)
            reloaded = ad.read_h5ad(path)
            assert_matrix(reloaded.X, result.X)
            equal(reloaded.obs, result.adata.obs)
            equal(reloaded.var, result.adata.var)
            for key in result.adata.layers:
                assert_matrix(reloaded.layers[key], result.adata.layers[key])
            equal(reloaded.obsm["X_embedding"], result.adata.obsm["X_embedding"])
            assert_matrix(reloaded.obsp["connectivities"], result.adata.obsp["connectivities"])
            assert reloaded.raw is None and not reloaded.varm and not reloaded.varp
            assert reloaded.uns["heart_orthomap"]["source_count_duplication"] == False
        for method in ("one2one", "sum", "best_sequence", "best_expression", "highest_total_expression"):
            self.case("T24_h5ad_roundtrip_" + method, lambda method=method: roundtrip(method))


def mutation_checks():
    golden = load_json("golden.json")["core"]
    caught = []
    mutations = {
        "Cartesian_count_inflation": [[6, 86, 0], [0, 12, 0], [12, 2, 0]],
        "sum_instead_of_best_sequence": golden["sum"],
        "per_cell_max_instead_of_global_choice": [[2, 40, 0], [0, 5, 0], [4, 1, 0]],
        "signed_int32_wrap": [[-2147483648]],
    }
    for name, values in mutations.items():
        try:
            if name == "signed_int32_wrap":
                assert_matrix(sparse.csr_matrix(values), [[2147483648]])
            else:
                method = "best_expression" if name.startswith("per_cell") else "best_sequence"
                assert_result(Result(sparse.csr_matrix(values), tuple(golden["human_ids"])), golden["human_ids"], golden[method])
        except AssertionError:
            caught.append(name)
    try:
        assert_result(Result(sparse.csr_matrix([[45], [6], [5]]), ("HUMAN_DUP",)), golden["human_ids"], golden["sum"])
    except AssertionError:
        caught.append("symbol_collapse")
    assert len(caught) == 5, caught
    return caught


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--self-check", action="store_true")
    mode.add_argument("--installed", action="store_true", help="Validate the installed heart_orthomap")
    parser.add_argument("--report-name", default=None)
    args = parser.parse_args()
    reports = ROOT / "reports"
    reports.mkdir(exist_ok=True)
    package = None
    hashes_before = None
    if args.self_check:
        mapper = independent_oracle.map_counts
        detected = mutation_checks()
        identity = {"mode": "independent_oracle_only"}
    else:
        import package_adapter as package
        hashes_before = package.source_hashes()
        mapper = package.map_counts
        detected = []
        identity = dict(package.PACKAGE_IDENTITY, version=package.api().__version__)
    suite = Suite(mapper, package)
    results = suite.run()
    hashes_after = package.source_hashes() if package else None
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "package_executed": not args.self_check,
        "identity": identity, "python": sys.executable,
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__},
        "passed": sum(item["status"] == "passed" for item in results),
        "failed": sum(item["status"] == "failed" for item in results),
        "results": results, "detected_harness_mutations": detected,
        "source_hashes_before": hashes_before, "source_hashes_after": hashes_after,
        "package_source_unchanged_during_run": hashes_before == hashes_after if package else None,
        "fixture_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT / "fixtures").glob("*.json"))},
        "unexecuted_release_checks": ["T21_all_type_filter_combinations_not_exposed_by_API", "T22_resource_isolated_stress_separate_report", "T23_complete_audit_schema"],
        "release_ready": False,
    }
    name = args.report_name or ("self_check.json" if args.self_check else "package_validation.json")
    if Path(name).name != name:
        raise ValueError("Report name must be a filename, not a path")
    path = reports / name
    path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print("SUMMARY", json.dumps({"passed": report["passed"], "failed": report["failed"], "report": str(path)}))
    return int(bool(report["failed"]) or (package is not None and hashes_before != hashes_after))


if __name__ == "__main__":
    raise SystemExit(main())
