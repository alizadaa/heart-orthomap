from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest
import warnings

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from heart_orthomap import convert_anndata, load_mapping, load_bundled_mapping
from heart_orthomap.cli import main


class FrozenRelease(unittest.TestCase):
    def test_exact_release_and_quarantines(self):
        mouse = load_bundled_mapping("mouse")
        fish = load_bundled_mapping("zebrafish")
        self.assertEqual(len(mouse.edges), 25768)
        self.assertEqual(len(fish.edges), 22361)
        self.assertEqual(len(fish.issues), 645)
        self.assertEqual(fish.issues.issue.value_counts().to_dict(), {
            "unresolved_or_multiple_stable_IDs": 638,
            "one_feature_ID_assigned_to_multiple_canonical_genes": 6,
            "invalid_Ensembl_identifier_namespace": 1})

    def test_mh6_has_curated_missing_score_and_is_retained(self):
        fish = load_bundled_mapping("zebrafish")
        edge = fish.edges.loc[fish.edges.target_symbol.eq("MYH6")]
        self.assertEqual(edge.source_id.tolist(), ["ENSDARG00000090637"])
        self.assertTrue(edge.sequence_score.isna().all())
        data = ad.AnnData(X=sparse.csr_matrix([[7], [11]]),
            obs=pd.DataFrame(index=["c1", "c2"]), var=pd.DataFrame(index=edge.source_id))
        result = convert_anndata(data, fish, counts_layer=None, strategy="best_sequence")
        self.assertEqual(result.adata.var_names.tolist(), ["ENSG00000197616"])
        np.testing.assert_array_equal(result.adata.X.toarray(), [[7], [11]])

    def test_cli_audit_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "audit"
            with redirect_stdout(io.StringIO()):
                main(["audit", "--species", "mouse", "--out", str(path)])
            self.assertEqual(json.loads((path / "report.json").read_text())["unique_edges"], 25768)
            with redirect_stdout(io.StringIO()), self.assertRaises(SystemExit):
                main(["audit", "--species", "mouse", "--out", str(path)])

    def test_pair_flags_not_shared_triplet_filter(self):
        frame = pd.DataFrame({
            "human_ensembl_id": ["ENSG00000000001"] * 2,
            "human_symbol": ["H1"] * 2,
            "mouse_ensembl_id": ["ENSMUSG00000000001"] * 2,
            "mouse_symbol": ["M1"] * 2,
            "human_to_mouse_final_mapping_selected": [True, True],
            "human_to_fish_final_mapping_selected": [False, False],
            "human_to_mouse_final_type": ["One2One"] * 2,
            "human_to_mouse_final_reverse_degree": [1] * 2,
            "human_to_mouse_final_forward_degree": [1] * 2,
            "human_mouse_Ensembl116_pct_id_target_to_query": [80] * 2,
            "human_mouse_Ensembl116_pct_id_query_to_target": [90] * 2})
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "triad.tsv"
            frame.to_csv(path, sep="\t", index=False)
            with warnings.catch_warnings(record=True) as caught:
                table = load_mapping(path, "mouse")
            self.assertTrue(caught)
            self.assertEqual(len(table.edges), 1)
            self.assertEqual(table.edges.sequence_score.iloc[0], 85)

    def test_overlapping_identifier_collisions_are_quarantined_atomically(self):
        rows = []
        for source, source_key, target, target_key in (
            ("ENSDARG99000000001", "ZFIN_A", "ENSG99000000001", "HUMAN_A"),
            ("ENSDARG99000000001", "ZFIN_B", "ENSG99000000002", "HUMAN_B"),
            ("ENSDARG99000000002", "ZFIN_C", "ENSG99000000001", "HUMAN_C"),
            ("ENSDARG99000000003", "ZFIN_SAFE", "ENSG99000000003", "HUMAN_SAFE"),
        ):
            rows.append({"fish_ensembl_ids": source, "fish_key": source_key,
                "fish_symbol": source_key, "left_ensembl_ids": target,
                "left_key": target_key, "left_symbol": target_key,
                "final_mapping_selected": True, "human_to_fish_final_type": "One2One",
                "final_reverse_degree": 1, "final_forward_degree": 1,
                "Ensembl116_pct_id_target_to_query": 80,
                "Ensembl116_pct_id_query_to_target": 80})
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "collisions.tsv"
            pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
            table = load_mapping(path, "zebrafish")
        self.assertEqual(table.edges.source_id.tolist(), ["ENSDARG99000000003"])
        self.assertEqual(len(table.issues), 3)


if __name__ == "__main__":
    unittest.main()
