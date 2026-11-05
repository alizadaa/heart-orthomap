#!/usr/bin/env python3
"""Build compact distributable pairwise projections; never modify source tables."""
import argparse
import json
from pathlib import Path

from heart_orthomap import load_mapping
from heart_orthomap.mapping import file_sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError("Use a new output directory; frozen snapshots are not overwritten")
    args.out.mkdir(parents=True)
    manifest = {"release": "2026-10-19", "policy": "HumanMouseFrozen_ZFINFishPriority_EnsemblFallback_SecondarySupport_v1",
        "code_license": "MIT", "data_terms": "Derived from Ensembl, ZFIN and supporting sources; upstream terms apply; see DATA_SOURCES.md",
        "not_a_unified_orthology_probability": True, "files": {}}
    for species, name in (("mouse", "human_mouse_final_mapping.tsv"), ("zebrafish", "human_fish_final_mapping.tsv")):
        table = load_mapping(args.release_dir / name, species)
        path = args.out / (species + "_to_human.tsv.gz")
        table.edges.to_csv(path, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
        issue_file = ""
        if len(table.issues):
            issue_file = species + "_unresolved_identifiers.tsv.gz"
            table.issues.to_csv(args.out / issue_file, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
        manifest["files"][species] = {"file": path.name, "sha256": file_sha256(path),
            "upstream_file": name, "upstream_sha256": table.provenance["mapping_sha256"],
            "unique_resolved_edges": len(table.edges), "unresolved_rows_excluded": len(table.issues),
            "unresolved_identifiers_file": issue_file,
            "sequence_score_definition": table.provenance["sequence_score_definition"],
            "sequence_version": "Ensembl116", "target_species": "human"}
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
