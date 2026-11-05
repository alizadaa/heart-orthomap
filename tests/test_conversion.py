import copy
import tempfile
from pathlib import Path
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from heart_orthomap import MappingTable, convert_anndata, expression_similarity_scores, load_mapping


def fixture():
    edges = pd.DataFrame([
        ("m1", "h1", "M1", "H1", 80), ("m2", "h1", "M2", "H1", 90),
        ("m3", "h2", "M3", "H2", 99), ("m4", "h3", "M4", "H3", 60),
        ("m4", "h4", "M4", "H4", 70)],
        columns=["source_id", "target_id", "source_symbol", "target_symbol", "sequence_score"])
    table = MappingTable(edges, "mouse")
    data = ad.AnnData(X=sparse.csr_matrix([[1, 2, 4, 8, 16], [2, 3, 5, 7, 11]]),
        obs=pd.DataFrame({"stage": ["E11.5", "E12.5"], "sample": ["a", "b"]}, index=["cell2", "cell1"]),
        var=pd.DataFrame(index=["m1", "m2", "m3", "m4", "unmapped"]))
    data.layers["counts"] = data.X.copy()
    data.layers["lognorm"] = data.X.astype(float).copy()
    data.obsm["X_umap"] = np.array([[0., 1.], [2., 3.]])
    return data, table


class ConversionTests(unittest.TestCase):
    def test_sum_sparse_no_source_duplication(self):
        data, table = fixture()
        before = data.copy()
        result = convert_anndata(data, table, strategy="sum", source_ambiguity="drop")
        np.testing.assert_array_equal(result.adata.X.toarray(), [[3, 4], [5, 5]])
        pd.testing.assert_frame_equal(result.adata.obs[before.obs.columns], before.obs)
        np.testing.assert_array_equal(result.adata.obsm["X_umap"], before.obsm["X_umap"])
        self.assertEqual((data.X != before.X).nnz, 0)
        self.assertEqual(result.report["aggregated_human_features"], 1)
        self.assertNotIn("lognorm", result.adata.layers)
        self.assertFalse(result.report["source_count_duplication"])

    def test_one2one_is_default(self):
        data, table = fixture()
        result = convert_anndata(data, table)
        self.assertEqual(result.adata.var_names.tolist(), ["h2"])
        np.testing.assert_array_equal(result.adata.X.toarray(), [[4], [5]])

    def test_ambiguity_requires_explicit_policy(self):
        data, table = fixture()
        with self.assertRaisesRegex(ValueError, "multiple human"):
            convert_anndata(data, table, strategy="sum")

    def test_best_sequence_selects_not_sums(self):
        data, table = fixture()
        result = convert_anndata(data, table, strategy="best_sequence", source_ambiguity="drop")
        np.testing.assert_array_equal(result.adata.X.toarray(), [[2, 4], [3, 5]])
        self.assertEqual(result.report["aggregated_human_features"], 0)

    def test_expression_is_independent_from_sequence(self):
        data, table = fixture()
        scores = pd.DataFrame({"source_id": ["m1", "m2"], "target_id": ["h1", "h1"],
                               "expression_similarity": [.9, .2]})
        result = convert_anndata(data, table, strategy="best_expression", source_ambiguity="drop", scores=scores)
        np.testing.assert_array_equal(result.adata.X.toarray(), [[1, 4], [2, 5]])
        self.assertIn("expression_scores_sha256", result.report)

    def test_numeric_string_scores_are_normalized_before_thresholding(self):
        data, table = fixture()
        scores = pd.DataFrame({"source_id": ["m1", "m2"], "target_id": ["h1", "h1"],
                               "expression_similarity": ["0.9", "0.2"]})
        before = scores.copy(deep=True)
        result = convert_anndata(data, table, strategy="best_expression",
            source_ambiguity="drop", scores=scores, min_score=.5)
        np.testing.assert_array_equal(result.adata.X.toarray(), [[1], [2]])
        pd.testing.assert_frame_equal(scores, before)

    def test_legacy_abundance_selection_not_similarity(self):
        data, table = fixture()
        result = convert_anndata(data, table, strategy="highest_total_expression", source_ambiguity="drop")
        np.testing.assert_array_equal(result.adata.X.toarray(), [[2, 4], [3, 5]])
        self.assertEqual(result.report["aggregated_human_features"], 0)

    def test_missing_or_tied_scores_not_arbitrarily_resolved(self):
        data, table = fixture()
        for value in (80, np.nan):
            frame = table.edges.copy()
            frame.loc[frame.source_id == "m2", "sequence_score"] = value
            changed = MappingTable(frame, "mouse")
            with self.assertRaises(ValueError):
                convert_anndata(data, changed, strategy="best_sequence", source_ambiguity="drop")
            result = convert_anndata(data, changed, strategy="best_sequence", source_ambiguity="drop", tie_policy="drop")
            self.assertEqual(result.adata.var_names.tolist(), ["h2"])

    def test_repeated_edges_do_not_double_counts(self):
        data, table = fixture()
        doubled = MappingTable(pd.concat([table.edges, table.edges]), "mouse")
        expected = convert_anndata(data, table, strategy="sum", source_ambiguity="drop")
        observed = convert_anndata(data, doubled, strategy="sum", source_ambiguity="drop")
        np.testing.assert_array_equal(expected.adata.X.toarray(), observed.adata.X.toarray())

    def test_conflicting_duplicate_edges_and_types_fail(self):
        _, table = fixture()
        different = table.edges.iloc[[0]].copy()
        different.sequence_score = 0
        with self.assertRaises(ValueError):
            MappingTable(pd.concat([table.edges, different]), "mouse")
        changed = table.edges.copy()
        changed.human_perspective_type = "One2One"
        with self.assertRaises(ValueError):
            MappingTable(changed, "mouse")

    def test_dense_and_h5ad_roundtrip(self):
        data, table = fixture()
        data.layers["counts"] = data.X.toarray()
        result = convert_anndata(data, table, strategy="sum", source_ambiguity="drop", output_ids="symbol")
        self.assertEqual(result.adata.var_names.tolist(), ["H1", "H2"])
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "mapped.h5ad"
            result.adata.write_h5ad(path)
            recovered = ad.read_h5ad(path)
            np.testing.assert_array_equal(recovered.X.toarray(), result.adata.X.toarray())
            pd.testing.assert_frame_equal(recovered.obs, result.adata.obs)

    def test_invalid_counts_and_ids_fail(self):
        data, table = fixture()
        for value in (-1., .5, np.nan, np.inf, float(2**53)):
            bad = data.copy()
            matrix = bad.X.astype(float).toarray()
            matrix[0, 0] = value
            bad.layers["counts"] = matrix
            with self.subTest(value=value), self.assertRaises(ValueError):
                convert_anndata(bad, table)
        bad = data.copy()
        bad.var_names = ["m1", "m1", "m3", "m4", "unmapped"]
        with self.assertRaises(ValueError):
            convert_anndata(bad, table)

    def test_symbol_and_feature_key_inputs(self):
        data, table = fixture()
        data.var["symbols"] = ["M1", "M2", "M3", "M4", "Unknown"]
        result = convert_anndata(data, table, feature_key="symbols", feature_type="symbol")
        self.assertEqual(result.adata.var_names.tolist(), ["h2"])

    def test_version_stripping_collision_fails(self):
        data, table = fixture()
        data.var_names = ["ENSMUSG1.1", "ENSMUSG1.2", "m3", "m4", "unmapped"]
        with self.assertRaises(ValueError):
            convert_anndata(data, table, strip_ensembl_versions=True)

    def test_constant_expression_profiles_missing_not_zero(self):
        source = ad.AnnData(X=sparse.csr_matrix([[1, 9], [1, 9], [1, 9]]),
            obs=pd.DataFrame(index=["a", "b", "c"]), var=pd.DataFrame(index=["m", "other"]))
        source.layers["counts"] = source.X.copy()
        reference = source.copy()
        reference.var_names = ["h", "other"]
        mapping = MappingTable(pd.DataFrame([("m", "h", "M", "H")], columns=["source_id", "target_id", "source_symbol", "target_symbol"]), "mouse")
        result = expression_similarity_scores(source, reference, mapping, source_groups=["g1", "g2", "g3"],
            reference_groups=["g1", "g2", "g3"], min_cells=1)
        self.assertTrue(np.isnan(result.expression_similarity.iloc[0]))

    def test_groups_and_score_provenance(self):
        data, table = fixture()
        with self.assertRaises(ValueError):
            expression_similarity_scores(data, data, table, source_groups=["g", "g"], reference_groups=["g", "g"])
        with self.assertRaises(ValueError):
            convert_anndata(data, table, strategy="best_expression", source_ambiguity="drop")


if __name__ == "__main__":
    unittest.main()
