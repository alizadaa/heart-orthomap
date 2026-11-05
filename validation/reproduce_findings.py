#!/usr/bin/env python3
"""Minimal standalone reproducers. Never edits or patches the package."""

import argparse
from pathlib import Path
import sys

sys.dont_write_bytecode = True
if not __debug__:
    raise RuntimeError("Assertions are required; run without -O or PYTHONOPTIMIZE")

import anndata as ad
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parent
from heart_orthomap import MappingTable, convert_anndata, load_mapping


def numeric_string_scores():
    data = ad.AnnData(sparse.csr_matrix([[1, 9]]), var=pd.DataFrame(index=["S_A", "S_B"]))
    mapping = MappingTable(pd.DataFrame({
        "source_id": ["S_A", "S_B"], "target_id": ["H", "H"],
        "source_symbol": ["A", "B"], "target_symbol": ["H", "H"],
    }), "mouse")
    scores = pd.DataFrame({"source_id": ["S_A", "S_B"], "target_id": ["H", "H"],
                           "expression_similarity": ["0.9", "0.1"]})
    result = convert_anndata(data, mapping, strategy="best_expression", counts_layer=None,
                             scores=scores, min_score=0.5)
    assert result.adata.X[0, 0] == 1
    print("PASS: numeric-string scores were normalized before thresholding")


def overlapping_frame():
    rows = []
    for source, source_key, target, target_key in (
        ("ENSDARG99000000001", "ZFIN_A", "ENSG99000000001", "HUMAN_A"),
        ("ENSDARG99000000001", "ZFIN_B", "ENSG99000000002", "HUMAN_B"),
        ("ENSDARG99000000002", "ZFIN_C", "ENSG99000000001", "HUMAN_C"),
        ("ENSDARG99000000003", "ZFIN_SAFE", "ENSG99000000003", "HUMAN_SAFE"),
    ):
        rows.append({"fish_ensembl_ids": source, "fish_key": source_key, "fish_symbol": source_key,
            "left_ensembl_ids": target, "left_key": target_key, "left_symbol": target_key,
            "final_mapping_selected": True, "human_to_fish_final_type": "One2One",
            "final_reverse_degree": 1, "final_forward_degree": 1,
            "Ensembl116_pct_id_target_to_query": 80, "Ensembl116_pct_id_query_to_target": 80})
    return pd.DataFrame(rows)


def overlapping_collisions():
    path = ROOT / "reports" / "overlapping_collision_fixture.txt"
    path.parent.mkdir(exist_ok=True)
    overlapping_frame().to_csv(path, sep="\t", index=False)
    mapping = load_mapping(path, "zebrafish")
    print("SURVIVING PAIRS", mapping.edges[["source_id", "target_id"]].to_dict("records"))
    print("QUARANTINED", mapping.issues[["source_id", "target_id", "issue"]].to_dict("records"))
    assert set(mapping.edges.source_id) == {"ENSDARG99000000003"}, "Y->H_SHARED survived despite H_SHARED identifying two canonical humans"
    assert len(mapping.issues) == 3
    print("PASS: collisions on both axes quarantined against the original selected graph")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("case", choices=("numeric-string-scores", "overlapping-collisions"))
    choice = parser.parse_args().case
    if choice == "numeric-string-scores":
        numeric_string_scores()
    else:
        overlapping_collisions()
