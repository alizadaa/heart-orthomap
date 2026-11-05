"""Run offline: python examples/convert_toy.py --out NEW_DIRECTORY."""
import argparse
import json
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy import sparse

from heart_orthomap import MappingTable, convert_anndata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError("Choose a new output directory")
    # Synthetic paralogs, not claims about real genes or their orthology.
    mapping = MappingTable(pd.DataFrame([
        ("mouse_a", "human_a", "MouseA", "HUMANA", 95.),
        ("mouse_b", "human_a", "MouseB", "HUMANA", 80.),
        ("mouse_c", "human_c", "MouseC", "HUMANC", 99.)],
        columns=["source_id", "target_id", "source_symbol", "target_symbol", "sequence_score"]), "mouse")
    data = ad.AnnData(X=sparse.csr_matrix([[2, 3, 7], [4, 1, 9]]),
        obs=pd.DataFrame({"stage": ["stage1", "stage2"]}, index=["barcode1", "barcode2"]),
        var=pd.DataFrame(index=["mouse_a", "mouse_b", "mouse_c"]))
    data.layers["counts"] = data.X.copy()
    args.out.mkdir(parents=True)
    summary = {}
    expression_scores = pd.DataFrame({"source_id": ["mouse_a", "mouse_b"],
        "target_id": ["human_a", "human_a"], "expression_similarity": [.3, .9]})
    for strategy in ("one2one", "sum", "best_sequence", "best_expression", "highest_total_expression"):
        result = convert_anndata(data, mapping, strategy=strategy,
            scores=expression_scores if strategy == "best_expression" else None)
        result.adata.write_h5ad(args.out / (strategy + ".h5ad"))
        result.decisions.to_csv(args.out / (strategy + "_decisions.tsv"), sep="\t", index=False)
        summary[strategy] = {"genes": result.adata.var_names.tolist(), "counts": result.adata.X.toarray().tolist()}
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
