"""Read frozen pairwise tables or project explicit decisions from a three-species table."""
from dataclasses import dataclass, field
import hashlib
from importlib.resources import files
from pathlib import Path
import warnings

import numpy as np
import pandas as pd

TYPES = {"One2One", "One2Many", "Many2One", "Many2Many"}
REQUIRED = ("source_id", "target_id", "source_symbol", "target_symbol")


@dataclass
class MappingTable:
    edges: pd.DataFrame
    source_species: str
    provenance: dict = field(default_factory=dict)
    issues: pd.DataFrame = field(default_factory=pd.DataFrame)

    def __post_init__(self):
        self.edges = _canonical(self.edges)


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1048576), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _selected(series):
    values = series.astype(str).str.lower()
    if not values.isin(["true", "false", "1", "0"]).all():
        raise ValueError("Selection flags must explicitly be True/False or 1/0")
    return values.isin(["true", "1"])


def _numeric(series):
    values = pd.to_numeric(series.replace("", np.nan), errors="raise").astype(float)
    if np.isinf(values).any():
        raise ValueError("Infinite mapping scores are invalid")
    return values


def _canonical(frame):
    frame = frame.copy(deep=True)
    if not set(REQUIRED).issubset(frame.columns):
        raise ValueError("Canonical mapping requires: " + ", ".join(REQUIRED))
    for key in REQUIRED:
        if frame[key].isna().any():
            raise ValueError("Missing mapping identifier/symbol: " + key)
        frame[key] = frame[key].astype(str)
    for key in ("source_id", "target_id"):
        if frame[key].str.contains(r"[;\t\r\n]", regex=True).any() or frame[key].eq("").any():
            raise ValueError("Each edge needs one resolved stable identifier: " + key)
    for key in ("sequence_score", "expression_similarity"):
        if key in frame:
            frame[key] = _numeric(frame[key])
    if "sequence_score" in frame:
        valid = frame.sequence_score.dropna()
        if ((valid < 0) | (valid > 100)).any():
            raise ValueError("Sequence scores must be percentage identities in [0,100]")
    if "expression_similarity" in frame:
        valid = frame.expression_similarity.dropna()
        if ((valid < -1) | (valid > 1)).any():
            raise ValueError("Expression-similarity scores must be correlations in [-1,1]")
    # Triad Cartesian repetitions never become repeated weights in a count projection.
    comparable = frame.mask(frame.eq(""))
    unique = comparable.groupby(["source_id", "target_id"], sort=False).nunique(dropna=True)
    if unique.gt(1).any().any():
        conflicts = unique.columns[unique.gt(1).any()].tolist()
        raise ValueError("Conflicting metadata for repeated stable-ID pair: " + ", ".join(conflicts))
    frame = comparable.groupby(["source_id", "target_id"], sort=True, as_index=False).first()
    frame[["source_symbol", "target_symbol"]] = frame[["source_symbol", "target_symbol"]].fillna("")
    if len(frame) == 0:
        raise ValueError("Mapping has no resolved selected pairs")
    observed_source = frame.groupby("source_id").target_id.transform("nunique").astype(int)
    observed_target = frame.groupby("target_id").source_id.transform("nunique").astype(int)
    for key, observed in (("source_degree", observed_source), ("target_degree", observed_target)):
        if key in frame:
            supplied = _numeric(frame[key])
            if supplied.isna().any() or (supplied < observed).any() or (supplied != np.floor(supplied)).any():
                raise ValueError("Supplied graph degrees cannot be below observed unique degrees")
            frame[key] = supplied.astype(int)
        else:
            frame[key] = observed
    expected_types = ["One2One" if s == t == 1 else "One2Many" if s == 1 else
                      "Many2One" if t == 1 else "Many2Many"
                      for s, t in zip(frame.source_degree, frame.target_degree)]
    if "human_perspective_type" not in frame:
        frame["human_perspective_type"] = expected_types
    if not frame.human_perspective_type.isin(TYPES).all():
        raise ValueError("Unknown final orthology type")
    if frame.human_perspective_type.tolist() != expected_types:
        raise ValueError("Final orthology types disagree with the complete-graph degrees")
    for id_key, symbol_key in (("source_id", "source_symbol"), ("target_id", "target_symbol")):
        if frame.groupby(id_key)[symbol_key].nunique().gt(1).any():
            raise ValueError("One stable ID has conflicting symbols: " + id_key)
    return frame.reset_index(drop=True)


def load_mapping(path, source_species="mouse", sequence_version="ensembl116"):
    """Read selected pairwise edges; the target is human, never an inferred mouse-fish link.

    Prefer the complete pairwise final table. Three-species tables are accepted
    with an explicit warning because their pairwise coverage can be incomplete.
    """
    species = {"fish": "zebrafish"}.get(source_species, source_species)
    if species not in {"mouse", "zebrafish"}:
        raise ValueError("Bundled readers currently support mouse/zebrafish to human")
    if sequence_version not in {"ensembl116", "archived"}:
        raise ValueError("sequence_version must be ensembl116 or archived")
    path = Path(path)
    raw = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    provenance = {"mapping_file": path.name, "mapping_sha256": file_sha256(path),
                  "source_species": species, "target_species": "human", "input_rows": len(raw),
                  "sequence_version": sequence_version}
    if set(REQUIRED).issubset(raw.columns):
        if "selected" in raw:
            raw = raw.loc[_selected(raw.selected)]
        provenance["input_schema"] = "canonical_pairwise"
        return MappingTable(raw, species, provenance)
    triad = "human_to_fish_final_mapping_selected" in raw
    if triad:
        warnings.warn("Three-species projection: use complete human_mouse_final_mapping.tsv or "
                      "human_fish_final_mapping.tsv for full pairwise coverage.", UserWarning, stacklevel=2)
        prefix = "human_mouse" if species == "mouse" else "human_fish"
        short = "mouse" if species == "mouse" else "fish"
        selected = _selected(raw[f"human_to_{short}_final_mapping_selected"])
        columns = dict(source_id="mouse_ensembl_id" if species == "mouse" else "fish_ensembl_ids",
            target_id="human_ensembl_id", source_symbol="mouse_symbol" if species == "mouse" else "fish_symbol",
            target_symbol="human_symbol", human_perspective_type=f"human_to_{short}_final_type",
            source_degree=f"human_to_{short}_final_reverse_degree", target_degree=f"human_to_{short}_final_forward_degree")
        if species == "zebrafish":
            columns["source_gene_key"] = "fish_key"
        score_prefix = prefix + ("_Ensembl116_" if sequence_version == "ensembl116" else "_Ensembl_")
    elif species == "mouse" and "human_ensembl_id" in raw:
        selected = _selected(raw.final_mapping_selected)
        columns = dict(source_id="mouse_ensembl_id", target_id="human_ensembl_id", source_symbol="mouse_symbol",
            target_symbol="human_symbol", human_perspective_type="human_to_mouse_final_type",
            source_degree="human_to_mouse_final_reverse_degree", target_degree="human_to_mouse_final_forward_degree")
        score_prefix = "Ensembl116_" if sequence_version == "ensembl116" else "archived_Ensembl_"
    elif species == "zebrafish" and "left_ensembl_ids" in raw:
        selected = _selected(raw.final_mapping_selected)
        columns = dict(source_id="fish_ensembl_ids", target_id="left_ensembl_ids", source_symbol="fish_symbol",
            target_symbol="left_symbol", human_perspective_type="human_to_fish_final_type",
            source_degree="final_reverse_degree", target_degree="final_forward_degree")
        columns.update(source_gene_key="fish_key", target_gene_key="left_key")
        score_prefix = "Ensembl116_"
        if sequence_version == "archived":
            raise ValueError("This pairwise fish table does not supply archived exact-pair scores")
    else:
        raise ValueError("Unrecognized final mapping schema/species")
    chosen = raw.loc[selected].copy()
    edges = pd.DataFrame({out: chosen[native] for out, native in columns.items()})
    usable = ~edges.source_id.eq("") & ~edges.target_id.eq("")
    usable &= ~edges.source_id.str.contains(";", regex=False) & ~edges.target_id.str.contains(";", regex=False)
    issue = pd.Series("", index=edges.index)
    issue.loc[~usable] = "unresolved_or_multiple_stable_IDs"
    source_pattern = r"ENSMUSG\d{11}" if species == "mouse" else r"ENSDARG\d{11}"
    invalid = usable & (~edges.source_id.str.fullmatch(source_pattern) | ~edges.target_id.str.fullmatch(r"ENSG\d{11}"))
    issue.loc[invalid] = "invalid_Ensembl_identifier_namespace"
    usable.loc[invalid] = False
    collision_scope = usable.copy()
    collisions = pd.Series(False, index=edges.index)
    for id_key, canonical_key in (("source_id", "source_gene_key"), ("target_id", "target_gene_key")):
        if canonical_key in edges:
            counts = edges.loc[collision_scope].groupby(id_key)[canonical_key].nunique()
            collided = collision_scope & edges[id_key].isin(counts.index[counts.gt(1)])
            issue.loc[collided] = "one_feature_ID_assigned_to_multiple_canonical_genes"
            collisions |= collided
    usable.loc[collisions] = False
    issues = edges.loc[~usable].assign(issue=issue.loc[~usable]).reset_index(drop=True)
    edges = edges.loc[usable].copy()
    score_names = [score_prefix + k for k in ("pct_id_target_to_query", "pct_id_query_to_target")]
    if not set(score_names).issubset(chosen.columns):
        raise ValueError("Requested sequence-score source is not present; no silent fallback")
    left, right = (_numeric(chosen.loc[usable, k]) for k in score_names)
    if ((left.dropna() < 0) | (left.dropna() > 100)).any() or ((right.dropna() < 0) | (right.dropna() > 100)).any():
        raise ValueError("Percentage identities outside [0,100]")
    edges["sequence_score"] = (left + right) / 2
    if "final_mapping_policy" in chosen:
        edges["mapping_policy"] = chosen.loc[usable, "final_mapping_policy"]
    reason_key = "human_to_fish_final_reason" if triad and species == "zebrafish" else "final_mapping_reason"
    if reason_key in chosen:
        edges["evidence"] = chosen.loc[usable, reason_key]
    provenance.update(input_schema="three_species_projection" if triad else "complete_pairwise_final",
                      selected_rows=len(chosen), excluded_unresolved_rows=len(issues),
                      sequence_score_definition="arithmetic mean of both exact-pair percentage identities; missing remains missing")
    return MappingTable(edges.reset_index(drop=True), species, provenance, issues)


def load_bundled_mapping(source_species="mouse"):
    species = {"fish": "zebrafish"}.get(source_species, source_species)
    if species not in {"mouse", "zebrafish"}:
        raise ValueError("Bundled species: mouse, zebrafish")
    root = files("heart_orthomap").joinpath("data")
    path = root.joinpath(species + "_to_human.tsv.gz")
    mapping = load_mapping(str(path), species)
    import json
    manifest = json.loads(root.joinpath("manifest.json").read_text())
    if mapping.provenance["mapping_sha256"] != manifest["files"][species]["sha256"]:
        raise ValueError("Bundled mapping checksum differs from its frozen release")
    mapping.provenance.update(manifest["files"][species])
    mapping.provenance["release"] = manifest["release"]
    issue_file = manifest["files"][species].get("unresolved_identifiers_file")
    if issue_file:
        mapping.issues = pd.read_csv(str(root.joinpath(issue_file)), sep="\t", dtype=str, keep_default_na=False)
    return mapping
