"""API-handoff extensions, still independent of the package implementation."""

import copy
import math
import warnings

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from fixture_support import ROOT, core_request, load_json, pair_edge, small_request
from validation_contract import Request


def expression_fixture():
    groups = np.repeat(["A", "B", "C"], 5)
    source_counts = np.asarray([
        [9, 1, 10], [1, 1, 98], [4, 3, 3], [0, 5, 5], [8, 2, 10],
        [1, 9, 10], [1, 10, 89], [2, 6, 2], [0, 8, 2], [1, 9, 10],
        [4, 4, 12], [3, 1, 96], [5, 2, 3], [0, 4, 6], [7, 1, 12],
    ], dtype=np.int64)
    reference_counts = np.asarray([
        [8, 12], [7, 93], [8, 2], [6, 4], [9, 11],
        [1, 19], [1, 99], [2, 8], [1, 9], [2, 18],
        [4, 16], [3, 97], [5, 5], [3, 7], [5, 15],
    ], dtype=np.int64)
    source_counts = np.column_stack([source_counts, np.zeros(15, dtype=np.int64)])
    reference_counts = np.column_stack([reference_counts, np.zeros(15, dtype=np.int64)])
    source = ad.AnnData(sparse.csr_matrix((15, 4)),
                       obs=pd.DataFrame(index=["src_%02d" % i for i in range(15)]),
                       var=pd.DataFrame(index=["S1", "S2", "S_BG", "S_CONST"]))
    reference = ad.AnnData(sparse.csr_matrix((15, 3)),
                          obs=pd.DataFrame(index=["ref_%02d" % i for i in range(15)]),
                          var=pd.DataFrame(index=["H1", "H_BG", "H_CONST"]))
    source.layers["counts"] = sparse.csr_matrix(source_counts)
    reference.layers["counts"] = sparse.csr_matrix(reference_counts)
    edges = small_request([[1, 1, 1]], ["S1", "S2", "S_CONST"], [
        pair_edge("H1", "S1", 90), pair_edge("H1", "S2", 80),
        pair_edge("H_CONST", "S_CONST", 90), pair_edge("H_ABSENT", "S_ABSENT", 90),
    ]).edges
    return source, reference, groups, edges


def scalar_profiles(counts, labels, feature):
    result = []
    for group in sorted(set(labels)):
        rows = np.flatnonzero(labels == group)
        normalized = [10000 * int(counts[row, feature]) / sum(int(x) for x in counts[row]) for row in rows]
        result.append(math.log1p(math.fsum(normalized) / len(rows)))
    return result


def scalar_pearson(x, y):
    left = [value - math.fsum(x) / len(x) for value in x]
    right = [value - math.fsum(y) / len(y) for value in y]
    denominator = math.sqrt(math.fsum(v*v for v in left) * math.fsum(v*v for v in right))
    return math.fsum(a*b for a, b in zip(left, right)) / denominator if denominator else float("nan")


def run_extended(suite, assert_result, assert_matrix, reject, equal):
    adapter = suite.package
    api = adapter.api()

    def highest_total():
        request = core_request("highest_total_expression")
        golden = load_json("golden.json")["core"]
        assert_result(suite.checked(request), golden["human_ids"], golden["best_expression"])
        request = small_request([[9, 1], [0, 7]], ["S_A", "S_B"], method="highest_total_expression")
        assert_result(suite.checked(request), ["H"], [[9], [0]])
    suite.case("T25_highest_total_global_after_ambiguity_drop", highest_total)

    def total_overflow():
        request = small_request([[2**52, 1], [2**52, 1]], method="highest_total_expression")
        reject(suite.mapper, request)
    suite.case("T25_highest_total_unsafe_total_rejected", total_overflow)

    def total_tie():
        request = small_request([[9, 1], [1, 9]], method="highest_total_expression")
        reject(suite.mapper, request)
    suite.case("T25_highest_total_tie_rejected", total_tie, "P1")

    def scores_numeric_strings():
        request = small_request([[1, 9]], method="best_expression")
        request.scores = [{"source_id": "S000", "target_id": "H", "expression_similarity": "0.9"},
                          {"source_id": "S001", "target_id": "H", "expression_similarity": "0.1"}]
        request.min_score = 0.5
        assert_result(suite.checked(request), ["H"], [[1]])
    suite.case("T26_numeric_string_scores_with_threshold", scores_numeric_strings, "P1")

    def scoring_matches_scalar():
        source, reference, groups, edges = expression_fixture()
        mapping = adapter.mapping_table(edges)
        source_before = source.layers["counts"].copy()
        reference_before = reference.layers["counts"].copy()
        result = api.expression_similarity_scores(source, reference, mapping,
            source_groups=pd.Series(groups, index=source.obs_names),
            reference_groups=pd.Series(groups, index=reference.obs_names))
        counts_left, counts_right = source_before.toarray(), reference_before.toarray()
        expected = [scalar_pearson(scalar_profiles(counts_left, groups, i),
                                  scalar_profiles(counts_right, groups, 0)) for i in (0, 1)]
        values = result.set_index(["source_id", "target_id"]).expression_similarity
        for source_id, value in zip(("S1", "S2"), expected):
            assert abs(values.loc[(source_id, "H1")] - value) < 1e-12
        assert math.isnan(values.loc[("S_CONST", "H_CONST")])
        assert len(result) == 3 and (result.matched_groups == 3).all()
        assert result.attrs["heart_orthomap_scoring"]["matched_groups"] == ["A", "B", "C"]
        assert_matrix(source.layers["counts"], source_before)
        assert_matrix(reference.layers["counts"], reference_before)
        assert source.X.nnz == reference.X.nnz == 0
        # Independent cell orders must not alter matched-group correlations.
        source_order = np.random.default_rng(91).permutation(15)
        reference_order = np.random.default_rng(92).permutation(15)
        permuted = api.expression_similarity_scores(source[source_order].copy(), reference[reference_order].copy(), mapping,
            source_groups=groups[source_order], reference_groups=groups[reference_order])
        np.testing.assert_allclose(result.expression_similarity, permuted.expression_similarity, rtol=0, atol=1e-12, equal_nan=True)
    suite.case("T27_expression_scalar_CP10K_mean_log_Pearson", scoring_matches_scalar)

    for attack in ("series_wrong_order", "array_wrong_length", "null_group", "insufficient_groups",
                   "zero_library", "min_cells_bool", "min_groups_two"):
        def invalid_scoring(attack=attack):
            source, reference, groups, edges = expression_fixture()
            kwargs = {"source_groups": groups.copy(), "reference_groups": groups.copy()}
            if attack == "series_wrong_order":
                kwargs["source_groups"] = pd.Series(groups, index=source.obs_names[::-1])
            elif attack == "array_wrong_length":
                kwargs["source_groups"] = groups[:-1]
            elif attack == "null_group":
                kwargs["source_groups"] = groups.astype(object)
                kwargs["source_groups"][0] = None
            elif attack == "insufficient_groups":
                kwargs["min_cells"] = 6
            elif attack == "zero_library":
                source.layers["counts"] = source.layers["counts"].tolil()
                source.layers["counts"][0, :] = 0
                source.layers["counts"] = source.layers["counts"].tocsr()
            elif attack == "min_cells_bool":
                kwargs["min_cells"] = True
            elif attack == "min_groups_two":
                kwargs["min_groups"] = 2
            try:
                api.expression_similarity_scores(source, reference, adapter.mapping_table(edges), **kwargs)
            except ValueError as error:
                assert str(error)
            else:
                raise AssertionError("Invalid scoring input silently accepted: " + attack)
        suite.case("T27_expression_reject_" + attack, invalid_scoring, "P1")

    def mouse_reader_frame():
        rows = []
        for source, left, right in (("ENSMUSG99000000001", 100, 0), ("ENSMUSG99000000002", 40, 40), ("ENSMUSG99000000003", "", 80)):
            rows.append({"human_ensembl_id": "ENSG99000000001", "mouse_ensembl_id": source,
                "human_symbol": "H_SYMBOL", "mouse_symbol": source, "final_mapping_selected": "True",
                "human_to_mouse_final_type": "One2Many", "human_to_mouse_final_reverse_degree": 1,
                "human_to_mouse_final_forward_degree": 3,
                "Ensembl116_pct_id_target_to_query": left, "Ensembl116_pct_id_query_to_target": right,
                "archived_Ensembl_pct_id_target_to_query": 99, "archived_Ensembl_pct_id_query_to_target": 99})
        return pd.DataFrame(rows)

    def reader_arithmetic_missing():
        frame = mouse_reader_frame()
        path = ROOT / "reports" / "reader_mouse_arithmetic.txt"
        frame.to_csv(path, sep="\t", index=False)
        mapping = api.load_mapping(path, "mouse")
        scores = mapping.edges.set_index("source_id").sequence_score
        assert scores["ENSMUSG99000000001"] == 50 and scores["ENSMUSG99000000002"] == 40
        assert pd.isna(scores["ENSMUSG99000000003"]), "Missing current score fell back to archived score"
        archived = api.load_mapping(path, "mouse", sequence_version="archived")
        assert (archived.edges.sequence_score == 99).all()
        request = small_request([[4, 7]], ["ENSMUSG99000000001", "ENSMUSG99000000002"])
        data = ad.AnnData(request.X, var=pd.DataFrame(index=request.source_ids))
        result = api.convert_anndata(data, mapping, strategy="best_sequence", counts_layer=None)
        assert_matrix(result.adata.X, [[4]])
    suite.case("T28_reader_arithmetic_mean_no_archived_fallback", reader_arithmetic_missing)

    def triplet_reader():
        frame = mouse_reader_frame().iloc[:2].copy()
        frame["human_to_mouse_final_mapping_selected"] = True
        frame["human_to_fish_final_mapping_selected"] = False
        frame["human_to_mouse_final_forward_degree"] = 2
        frame["human_mouse_Ensembl116_pct_id_target_to_query"] = frame.Ensembl116_pct_id_target_to_query
        frame["human_mouse_Ensembl116_pct_id_query_to_target"] = frame.Ensembl116_pct_id_query_to_target
        repeated = pd.concat([frame.assign(fish_ensembl_ids="FISH_%d" % index) for index in range(4)], ignore_index=True)
        path = ROOT / "reports" / "triplet_Cartesian_fixture.txt"
        repeated.to_csv(path, sep="\t", index=False)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            mapping = api.load_mapping(path, "mouse")
        assert len(mapping.edges) == 2
        assert any("Three-species" in str(w.message) for w in caught)
        data = ad.AnnData(sparse.csr_matrix([[4, 7]]), var=pd.DataFrame(index=["ENSMUSG99000000001", "ENSMUSG99000000002"]))
        result = api.convert_anndata(data, mapping, strategy="sum", counts_layer=None)
        assert_matrix(result.adata.X, [[11]])
    suite.case("T02_native_triplet_Cartesian_reader", triplet_reader)

    def fish_collisions():
        rows = []
        collisions = {"ENSDARG99000000101": ("SYN_A1", "SYN_A2"),
                      "ENSDARG99000000102": ("SYN_B1", "SYN_B2"),
                      "ENSDARG99000000103": ("SYN_C1", "SYN_C2")}
        human_number = 0
        for source, symbols in dict(collisions, ENSDARG99000000001=("SAFE",)).items():
            for symbol in symbols:
                human_number += 1
                rows.append({"fish_ensembl_ids": source, "fish_key": "ZFIN_" + symbol, "fish_symbol": symbol,
                    "left_ensembl_ids": "ENSG%011d" % (99000000000 + human_number), "left_key": "HKEY_" + symbol, "left_symbol": symbol.upper(),
                    "final_mapping_selected": True, "human_to_fish_final_type": "One2One",
                    "final_reverse_degree": 1, "final_forward_degree": 1,
                    "Ensembl116_pct_id_target_to_query": 80, "Ensembl116_pct_id_query_to_target": 60})
        path = ROOT / "reports" / "fish_collision_fixture.txt"
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
        mapping = api.load_mapping(path, "zebrafish")
        assert set(mapping.edges.source_id) == {"ENSDARG99000000001"}
        assert len(mapping.issues) == 6 and set(mapping.issues.source_id) == set(collisions)
        assert set(mapping.issues.issue) == {"one_feature_ID_assigned_to_multiple_canonical_genes"}
        assert mapping.provenance["excluded_unresolved_rows"] == 6
    suite.case("T29_fish_canonical_ID_collision_quarantine", fish_collisions)

    def overlapping_collisions():
        from reproduce_findings import overlapping_frame
        path = ROOT / "reports" / "overlapping_collision_fixture.txt"
        overlapping_frame().to_csv(path, sep="\t", index=False)
        mapping = api.load_mapping(path, "zebrafish")
        assert set(mapping.edges.source_id) == {"ENSDARG99000000003"}, "Target ID collision disappeared after source-side quarantine"
        assert len(mapping.issues) == 3
    suite.case("T29_overlapping_source_target_collision_quarantine", overlapping_collisions, "P1")

    def bundled_species():
        collisions = {"ENSDARG00000041570", "ENSDARG00000038976", "ENSDARG00000087657"}
        summary = {}
        for species in ("mouse", "zebrafish"):
            mapping = api.load_bundled_mapping(species)
            assert not mapping.edges.duplicated(["source_id", "target_id"]).any()
            if species == "zebrafish":
                assert not (set(mapping.edges.source_id) & collisions)
                assert collisions.issubset(set(mapping.issues.source_id))
                assert len(mapping.issues.loc[mapping.issues.source_id.isin(collisions)]) == 6
            eligible = mapping.edges.loc[mapping.edges.source_degree.eq(1)].head(10)
            source_ids = sorted(eligible.source_id.unique())
            data = ad.AnnData(sparse.csr_matrix([list(range(1, len(source_ids)+1))]),
                             var=pd.DataFrame(index=source_ids))
            result = api.convert_anndata(data, mapping, strategy="sum", counts_layer=None)
            assert int(result.adata.X.sum()) == sum(range(1, len(source_ids)+1))
            summary[species] = {"edges": len(mapping.edges), "issues": len(mapping.issues),
                                "mapping_sha256": mapping.provenance["mapping_sha256"]}
        (ROOT / "reports" / "bundled_species_readback.json").write_text(__import__("json").dumps(summary, indent=2) + "\n")
    suite.case("T24_T29_bundled_mouse_fish_counts_and_collision_audit", bundled_species)
