"""Command-line count conversion and frozen-table audit."""
import argparse
import json
from pathlib import Path

import anndata as ad
import pandas as pd

from .convert import convert_anndata
from .mapping import load_mapping, load_bundled_mapping


def main(argv=None):
    parser = argparse.ArgumentParser(prog="heart-orthomap")
    parser.add_argument("--version", action="version", version="heart-orthomap 0.1.0")
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("convert", "audit"):
        command = sub.add_parser(name)
        command.add_argument("--species", choices=("mouse", "zebrafish"), required=True)
        command.add_argument("--mapping", type=Path, help="Final TSV/TSV.GZ; default: frozen bundled pairwise release")
        command.add_argument("--out", type=Path, required=True, help="New output directory, never overwritten")
    convert = sub.choices["convert"]
    convert.add_argument("--input", type=Path, required=True)
    convert.add_argument("--strategy", choices=("one2one", "sum", "best_sequence", "best_expression", "highest_total_expression"), default="one2one")
    convert.add_argument("--counts-layer", default="counts", help="Raw layer; use X explicitly for raw .X")
    convert.add_argument("--feature-key")
    convert.add_argument("--feature-type", choices=("ensembl", "symbol"), default="ensembl")
    convert.add_argument("--output-ids", choices=("ensembl", "symbol"), default="ensembl")
    convert.add_argument("--source-ambiguity", choices=("error", "drop"), default="error")
    convert.add_argument("--tie-policy", choices=("error", "drop"), default="error")
    convert.add_argument("--scores", type=Path, help="Matched-group expression scores TSV")
    convert.add_argument("--min-score", type=float)
    convert.add_argument("--strip-ensembl-versions", action="store_true")
    args = parser.parse_args(argv)
    if args.out.exists():
        parser.error("Output directory exists; choose a new path")
    mapping = load_mapping(args.mapping, args.species) if args.mapping else load_bundled_mapping(args.species)
    if args.command == "audit":
        args.out.mkdir(parents=True)
        mapping.edges.to_csv(args.out / "resolved_selected_edges.tsv.gz", sep="\t", index=False)
        mapping.issues.to_csv(args.out / "unresolved_identifiers.tsv", sep="\t", index=False)
        report = {**mapping.provenance, "unique_edges": len(mapping.edges),
            "human_genes": mapping.edges.target_id.nunique(), "source_genes": mapping.edges.source_id.nunique(),
            "edge_counts_by_human_perspective_type": mapping.edges.human_perspective_type.value_counts().to_dict(),
            "human_targets_with_multiple_sources": int(mapping.edges.groupby("target_id").source_id.nunique().gt(1).sum()),
            "sources_with_multiple_human_targets": int(mapping.edges.groupby("source_id").target_id.nunique().gt(1).sum())}
    else:
        data = ad.read_h5ad(args.input)
        result = convert_anndata(data, mapping, strategy=args.strategy,
            counts_layer=None if args.counts_layer == "X" else args.counts_layer,
            feature_key=args.feature_key, feature_type=args.feature_type, output_ids=args.output_ids,
            source_ambiguity=args.source_ambiguity, tie_policy=args.tie_policy,
            scores=pd.read_csv(args.scores, sep="\t") if args.scores else None,
            min_score=args.min_score, strip_ensembl_versions=args.strip_ensembl_versions)
        args.out.mkdir(parents=True)
        result.adata.write_h5ad(args.out / "human_mapped.h5ad", compression="gzip")
        result.decisions.to_csv(args.out / "mapping_decisions.tsv.gz", sep="\t", index=False)
        report = result.report
    (args.out / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": "completed", "output": str(args.out)}))


if __name__ == "__main__":
    main()
