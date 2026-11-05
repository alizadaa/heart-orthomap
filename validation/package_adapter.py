"""Thin schema translation to the installed package. No count behavior is repaired."""

import hashlib
from importlib import import_module
from importlib.resources import files

from validation_contract import Result

ADAPTER_READY = True
PACKAGE_IDENTITY = {"distribution": "heart-orthomap", "binding": "installed_public_api"}


def source_hashes():
    hashes = {}

    def walk(directory, prefix):
        for item in sorted(directory.iterdir(), key=lambda entry: entry.name):
            name = prefix + "/" + item.name
            if item.is_dir():
                walk(item, name)
            elif item.name.endswith(".py"):
                hashes[name] = hashlib.sha256(item.read_bytes()).hexdigest()

    walk(files(api()), "heart_orthomap")
    return hashes


def api():
    return import_module("heart_orthomap")


def mapping_table(edges):
    import pandas as pd

    rows = []
    names = {"one2one": "One2One", "one2many": "One2Many",
             "many2one": "Many2One", "many2many": "Many2Many"}
    for edge in edges:
        row = {
            "source_id": edge["source_id"], "target_id": edge["human_id"],
            "source_symbol": edge.get("source_symbol", edge["source_id"]),
            "target_symbol": edge.get("human_symbol", edge["human_id"]),
            "sequence_score": edge.get("sequence_score"),
            "human_perspective_type": names[edge.get("orthology_type", "one2many").lower()],
        }
        for key in ("source_degree", "target_degree"):
            if key in edge:
                row[key] = edge[key]
        rows.append(row)
    columns = ["source_id", "target_id", "source_symbol", "target_symbol", "sequence_score", "human_perspective_type"]
    frame = pd.DataFrame(rows) if rows else pd.DataFrame(columns=columns)
    return api().MappingTable(frame, "mouse", {"fixture": "independent_synthetic"})


def map_counts(request):
    import anndata as ad
    import pandas as pd

    data = request.adata
    if data is None:
        data = ad.AnnData(X=request.X, var=pd.DataFrame(index=list(request.source_ids)))
    scores = None if request.scores is None else pd.DataFrame(request.scores)
    result = api().convert_anndata(
        data, mapping_table(request.edges), strategy=request.method,
        counts_layer=request.count_layer, feature_key=request.feature_key,
        feature_type=request.feature_type, output_ids=request.output_ids,
        source_ambiguity=request.ambiguity, tie_policy=request.tie_policy,
        scores=scores, min_score=request.min_score,
        strip_ensembl_versions=request.strip_ensembl_versions,
        additional_count_layers=request.additional_count_layers,
    )
    selected = result.decisions.loc[result.decisions.selected]
    sources = {target: tuple(sorted(group.source_id))
               for target, group in selected.groupby("target_id", sort=True)}
    return Result(
        X=result.adata.X, human_ids=tuple(result.adata.var.human_ensembl_id),
        adata=result.adata, selected_sources=sources,
        provenance=result.report, decisions=result.decisions,
    )
