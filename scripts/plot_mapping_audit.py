#!/usr/bin/env python3
"""Reproduce descriptive figures from the frozen selected pairwise mappings.

No triplet table is read. Canonical fish keys are counted without expanding
semicolon-delimited feature-ID fields. Run with --self-test for small invariants.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "heart-orthomap-mpl"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RELEASE = (
    ROOT.parent
    / "outputs/orthology_human_mouse_zebrafish_ZFIN_2026_10_19"
    / "final_mapping_2026_10_19/release"
)
CARDINALITIES = ("One2One", "One2Many", "Many2One", "Many2Many")
CARD_LABELS = ("1:1", "1:many", "many:1", "many:many")
CARD_COLORS = ("#4477AA", "#CCBB44", "#AA3377", "#228833")
SYMBOLS = ("same_after_uppercase", "different_after_uppercase", "missing_either_symbol")
SYMBOL_LABELS = ("Same after uppercase", "Different after uppercase", "Missing either symbol")
SYMBOL_COLORS = ("#228877", "#CC6677", "#B8BCC2")
EVIDENCE = ("ZFIN_curated_primary", "Ensembl_fallback_no_ZFIN_ortholog_for_human_gene")
EVIDENCE_LABELS = ("ZFIN-curated primary", "Ensembl fallback")
EVIDENCE_COLORS = ("#4477AA", "#DD9944")
ID_STATUSES = ("single", "missing", "multiple", "invalid_syntax")
ID_LABELS = ("One valid-format ID", "No Ensembl ID", "Multiple Ensembl IDs", "Invalid ID syntax")
ID_COLORS = ("#228877", "#B8BCC2", "#AA3377", "#CC6677")
INK = "#242A30"
GREY = "#5F6871"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def boolean(value: str) -> bool:
    require(value in ("True", "False"), f"Unexpected boolean: {value!r}")
    return value == "True"


def feature_ids(value: str) -> frozenset[str]:
    """Inspect field coverage only; never turn these tokens into graph edges."""
    return frozenset(token.strip() for token in value.split(";") if token.strip())


def id_status(value: str, prefix: str) -> str:
    ids = feature_ids(value)
    if any(not re.fullmatch(prefix + r"\d+", item) for item in ids):
        return "invalid_syntax"
    n = len(ids)
    return "missing" if n == 0 else "single" if n == 1 else "multiple"


def symbol_status(human: str, source: str) -> str:
    human, source = human.strip(), source.strip()
    if not human or not source:
        return SYMBOLS[2]
    return SYMBOLS[0] if human.upper() == source.upper() else SYMBOLS[1]


def cardinality(forward: int, reverse: int) -> str:
    require(forward >= 1 and reverse >= 1, "An edge must have positive endpoint degrees")
    if forward == 1:
        return "One2One" if reverse == 1 else "Many2One"
    return "One2Many" if reverse == 1 else "Many2Many"


@dataclass
class Graph:
    name: str
    title: str
    path: Path
    raw_rows: int
    rows: list[dict[str, str]]
    human_column: str
    source_column: str
    human_symbol_column: str
    source_symbol_column: str
    cardinality_column: str
    human_neighbors: dict[str, set[str]]
    source_neighbors: dict[str, set[str]]

    @property
    def edge_count(self) -> int:
        return len(self.rows)

    def human_keys(self, rows=None) -> set[str]:
        return {r[self.human_column] for r in (self.rows if rows is None else rows)}

    def source_keys(self, rows=None) -> set[str]:
        return {r[self.source_column] for r in (self.rows if rows is None else rows)}

    def spelling(self, row: dict[str, str]) -> str:
        return symbol_status(row[self.human_symbol_column], row[self.source_symbol_column])


def read_graph(path: Path, species: str) -> Graph:
    if species == "mouse":
        h, s, hs, ss, c, fd, rd = (
            "human_ensembl_id", "mouse_ensembl_id", "human_symbol", "mouse_symbol",
            "human_to_mouse_final_type", "human_to_mouse_final_forward_degree",
            "human_to_mouse_final_reverse_degree",
        )
        extra = {"final_mapping_reason", "Ensembl116_exact_pair_status"}
    else:
        h, s, hs, ss, c, fd, rd = (
            "left_key", "fish_key", "left_symbol", "fish_symbol", "final_mapping_type",
            "final_forward_degree", "final_reverse_degree",
        )
        extra = {
            "left_ensembl_ids", "fish_ensembl_ids", "final_mapping_reason",
            "ZFIN_supported", "Ensembl_supported", "usable_for_Ensembl_feature_mapping",
            "requires_feature_ID_review", "Ensembl116_exact_pair_status",
            "human_to_fish_final_type",
        }
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {h, s, hs, ss, c, fd, rd, "final_mapping_selected"} | extra
        require(required <= set(reader.fieldnames or []), f"Missing columns in {path.name}: {required - set(reader.fieldnames or [])}")
        raw = list(reader)
    require(bool(raw), f"Empty mapping: {path}")
    edges = {}
    for row in raw:
        require(None not in row and all(v is not None for v in row.values()), "Malformed TSV row")
        require(boolean(row["final_mapping_selected"]), "Input contains a nonselected edge")
        key = (row[h], row[s])
        require(all(key), "Missing canonical endpoint key")
        require(key not in edges or edges[key] == row, f"Conflicting duplicate edge: {key}")
        edges[key] = row
    rows = [edges[key] for key in sorted(edges)]
    hn, sn = defaultdict(set), defaultdict(set)
    for row in rows:
        hn[row[h]].add(row[s])
        sn[row[s]].add(row[h])
    for row in rows:
        forward, reverse = len(hn[row[h]]), len(sn[row[s]])
        require((int(row[fd]), int(row[rd])) == (forward, reverse), f"Stored degrees disagree: {row[h]}, {row[s]}")
        require(row[c] == cardinality(forward, reverse), f"Stored cardinality disagrees: {row[h]}, {row[s]}")
        if species == "mouse":
            require(bool(re.fullmatch(r"ENSG\d+", row[h])), "Invalid human stable ID")
            require(bool(re.fullmatch(r"ENSMUSG\d+", row[s])), "Invalid mouse stable ID")
        else:
            require(row["human_to_fish_final_type"] == row[c], "Fish cardinality columns disagree")
            require(row["final_mapping_reason"] in EVIDENCE, "Unknown fish selection reason")
            curated = row["final_mapping_reason"] == EVIDENCE[0]
            require(boolean(row["ZFIN_supported"]) == curated, "ZFIN support conflicts with selection reason")
            require(curated or boolean(row["Ensembl_supported"]), "Fallback edge lacks archived Ensembl support")
            left = feature_ids(row["left_ensembl_ids"])
            fish = feature_ids(row["fish_ensembl_ids"])
            single_token = len(left) == 1 and len(fish) == 1
            require(boolean(row["usable_for_Ensembl_feature_mapping"]) == single_token, "Recorded usability disagrees with raw ID-token counts")
            require(boolean(row["requires_feature_ID_review"]) == (not single_token), "Recorded review flag disagrees with raw ID-token counts")
    return Graph(species, "Mouse" if species == "mouse" else "Zebrafish", path, len(raw), rows, h, s, hs, ss, c, hn, sn)


def unique_field_by_key(graph: Graph, key: str, field: str) -> dict[str, str]:
    values = defaultdict(set)
    for row in graph.rows:
        values[row[key]].add(row[field])
    require(all(len(v) == 1 for v in values.values()), f"Inconsistent {field} within canonical keys")
    return {k: next(iter(v)) for k, v in sorted(values.items())}


def mixed_cardinality_nodes(graph: Graph, column: str) -> int:
    memberships = defaultdict(set)
    for row in graph.rows:
        memberships[row[column]].add(row[graph.cardinality_column])
    return sum(len(value) > 1 for value in memberships.values())


def fraction(n: int, denominator: int) -> float:
    return n / denominator if denominator else 0.0


def build_tables(graphs: list[Graph]) -> dict[str, list[dict]]:
    tables = {name: [] for name in (
        "graph_totals", "cardinality_counts", "symbol_edge_counts", "symbol_human_key_counts",
        "fish_evidence_counts", "fish_id_coverage", "fish_feature_mapping_counts", "myh6_edges",
        "id_validation_issues", "fish_joint_id_coverage",
    )}
    for graph in graphs:
        tables["graph_totals"].append({
            "species_pair": f"human_{graph.name}", "input_rows": graph.raw_rows,
            "unique_edges": graph.edge_count, "exact_duplicate_rows_removed": graph.raw_rows - graph.edge_count,
            "unique_human_gene_keys": len(graph.human_neighbors),
            "unique_source_gene_keys": len(graph.source_neighbors),
            "human_keys_in_multiple_cardinality_classes": mixed_cardinality_nodes(graph, graph.human_column),
            "source_keys_in_multiple_cardinality_classes": mixed_cardinality_nodes(graph, graph.source_column),
            "human_keys_without_single_valid_Ensembl_ID": 0 if graph.name == "mouse" else sum(
                id_status(v, "ENSG") != "single" for v in unique_field_by_key(graph, graph.human_column, "left_ensembl_ids").values()
            ),
            "human_symbol_only_placeholder_keys": sum(key.startswith("human:SYMBOL:") for key in graph.human_neighbors),
        })
        for card in ("All",) + CARDINALITIES:
            rows = graph.rows if card == "All" else [r for r in graph.rows if r[graph.cardinality_column] == card]
            tables["cardinality_counts"].append({
                "species_pair": f"human_{graph.name}", "final_cardinality": card,
                "unique_edges": len(rows), "unique_human_gene_keys": len(graph.human_keys(rows)),
                "unique_source_gene_keys": len(graph.source_keys(rows)),
                "gene_counts_additive_across_cardinalities": False,
            })
            symbols = Counter(graph.spelling(row) for row in rows)
            for status in SYMBOLS:
                tables["symbol_edge_counts"].append({
                    "species_pair": f"human_{graph.name}", "final_cardinality": card,
                    "symbol_comparison": status, "unique_edges": symbols[status],
                    "denominator_unique_edges_in_class": len(rows),
                    "fraction_of_class_edges": fraction(symbols[status], len(rows)),
                    "denominator_unique_edges_in_graph": graph.edge_count,
                    "fraction_of_graph_edges": fraction(symbols[status], graph.edge_count),
                })
        human_symbols = defaultdict(set)
        for row in graph.rows:
            human_symbols[row[graph.human_column]].add(graph.spelling(row))
        human_counts = Counter()
        for states in human_symbols.values():
            if SYMBOLS[2] in states:
                label = "missing_symbol_on_at_least_one_edge"
            elif len(states) > 1:
                label = "mixed_same_and_different"
            elif states == {SYMBOLS[0]}:
                label = "all_edges_same_after_uppercase"
            else:
                label = "all_edges_different_after_uppercase"
            human_counts[label] += 1
        for label in ("all_edges_same_after_uppercase", "all_edges_different_after_uppercase", "mixed_same_and_different", "missing_symbol_on_at_least_one_edge"):
            tables["symbol_human_key_counts"].append({
                "species_pair": f"human_{graph.name}", "human_key_symbol_category": label,
                "unique_human_gene_keys": human_counts[label],
                "denominator_unique_human_gene_keys": len(human_symbols),
                "fraction_of_human_keys": fraction(human_counts[label], len(human_symbols)),
            })
        for row in graph.rows:
            # MYH6 is located by its stable human ID, not by a symbol join.
            stable = row[graph.human_column] if graph.name == "mouse" else row["left_ensembl_ids"]
            if stable != "ENSG00000197616":
                continue
            tables["myh6_edges"].append({
                "species_pair": f"human_{graph.name}", "human_key": row[graph.human_column],
                "human_ensembl_ids": stable, "human_symbol": row[graph.human_symbol_column],
                "source_key": row[graph.source_column], "source_symbol": row[graph.source_symbol_column],
                "source_ensembl_ids": row[graph.source_column] if graph.name == "mouse" else row["fish_ensembl_ids"],
                "final_cardinality": row[graph.cardinality_column],
                "selection_reason": row["final_mapping_reason"],
                "ZFIN_supported": row.get("ZFIN_supported", "not_applicable"),
                "archived_Ensembl_supported": row.get("Ensembl_supported", "not_applicable"),
                "Ensembl116_exact_pair_status": row["Ensembl116_exact_pair_status"],
                "symbol_comparison": graph.spelling(row),
            })
    fish = next(g for g in graphs if g.name == "fish")
    human_reasons = unique_field_by_key(fish, fish.human_column, "final_mapping_reason")
    for card in ("All",) + CARDINALITIES:
        rows = fish.rows if card == "All" else [r for r in fish.rows if r[fish.cardinality_column] == card]
        for reason in EVIDENCE:
            subset = [r for r in rows if r["final_mapping_reason"] == reason]
            tables["fish_evidence_counts"].append({
                "final_cardinality": card, "selection_reason": reason,
                "unique_edges": len(subset), "denominator_unique_edges_in_class": len(rows),
                "unique_human_gene_keys": len(fish.human_keys(subset)),
                "unique_source_gene_keys": len(fish.source_keys(subset)),
                "fraction_of_class_edges": fraction(len(subset), len(rows)),
            })
    require(sum(Counter(human_reasons.values()).values()) == len(fish.human_neighbors), "Human evidence partition failed")
    for side, field, key, prefix in (
        ("human", "left_ensembl_ids", fish.human_column, "ENSG"),
        ("fish", "fish_ensembl_ids", fish.source_column, "ENSDARG"),
    ):
        values = unique_field_by_key(fish, key, field)
        for unit, statuses in (
            ("unique_edges", Counter(id_status(r[field], prefix) for r in fish.rows)),
            ("unique_gene_keys", Counter(id_status(v, prefix) for v in values.values())),
        ):
            for status in ID_STATUSES:
                tables["fish_id_coverage"].append({
                    "side": side, "counting_unit": unit, "Ensembl_ID_status": status,
                    "count": statuses[status], "denominator": sum(statuses.values()),
                    "fraction": fraction(statuses[status], sum(statuses.values())),
                })
    for row in fish.rows:
        human_status = id_status(row["left_ensembl_ids"], "ENSG")
        fish_status = id_status(row["fish_ensembl_ids"], "ENSDARG")
        if "invalid_syntax" in (human_status, fish_status):
            tables["id_validation_issues"].append({
                "human_key": row[fish.human_column], "human_symbol": row[fish.human_symbol_column],
                "fish_key": row[fish.source_column], "fish_symbol": row[fish.source_symbol_column],
                "left_ensembl_ids": row["left_ensembl_ids"], "fish_ensembl_ids": row["fish_ensembl_ids"],
                "human_ID_status": human_status, "fish_ID_status": fish_status,
                "recorded_usable_for_Ensembl_feature_mapping": boolean(row["usable_for_Ensembl_feature_mapping"]),
                "syntax_audited_usable": human_status == fish_status == "single",
                "issue": "Non-Ensembl identifier in an Ensembl feature-ID field; input preserved unchanged",
            })
    joint = Counter((id_status(r["left_ensembl_ids"], "ENSG"), id_status(r["fish_ensembl_ids"], "ENSDARG")) for r in fish.rows)
    for human_status in ID_STATUSES:
        for fish_status in ID_STATUSES:
            tables["fish_joint_id_coverage"].append({
                "human_Ensembl_ID_status": human_status, "fish_Ensembl_ID_status": fish_status,
                "unique_edges": joint[human_status, fish_status], "denominator_unique_edges": fish.edge_count,
            })
    for reason in ("All",) + EVIDENCE:
        rows = fish.rows if reason == "All" else [r for r in fish.rows if r["final_mapping_reason"] == reason]
        for definition in ("recorded_flags", "syntax_audited"):
            for usable in (True, False):
                flags = [boolean(r["usable_for_Ensembl_feature_mapping"]) if definition == "recorded_flags" else (
                    id_status(r["left_ensembl_ids"], "ENSG") == id_status(r["fish_ensembl_ids"], "ENSDARG") == "single"
                ) for r in rows]
                count = sum(flag == usable for flag in flags)
                tables["fish_feature_mapping_counts"].append({
                    "selection_reason": reason, "coverage_definition": definition,
                    "usable_for_Ensembl_feature_mapping": usable,
                    "unique_edges": count, "denominator_unique_edges": len(rows),
                    "fraction_of_edges": fraction(count, len(rows)),
                })
    return tables


def style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 12,
        "axes.labelsize": 10, "axes.titleweight": "semibold", "text.color": INK,
        "axes.labelcolor": INK, "xtick.color": GREY, "ytick.color": INK,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.spines.left": False, "axes.edgecolor": "#CDD2D6",
        "axes.axisbelow": True, "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none", "savefig.facecolor": "white",
    })


def heading(fig, title: str, subtitle: str) -> None:
    fig.text(0.055, 0.97, title, fontsize=16, weight="bold", va="top")
    fig.text(0.055, 0.97 - 0.34 / fig.get_figheight(), subtitle, fontsize=10, color=GREY, va="top")


def footer(fig, text: str) -> None:
    fig.text(0.055, 0.025, text, fontsize=8.5, color=GREY, va="bottom", linespacing=1.5)


def save(fig, out: Path, name: str) -> None:
    # Stable PDF metadata avoids nondeterministic wall-clock timestamps.
    fig.savefig(out / f"{name}.png", dpi=300)
    fig.savefig(out / f"{name}.pdf", metadata={
        "Title": name, "Author": "Riley lab / Azad Alizada / CellTalk HHD",
        "Creator": "heart-orthomap plot_mapping_audit.py",
        "CreationDate": None, "ModDate": None,
    })
    plt.close(fig)


def legend_patches(labels, colors):
    return [Patch(facecolor=color, label=label) for label, color in zip(labels, colors)]


def plot_symbols(graphs: list[Graph], out: Path) -> None:
    fig = plt.figure(figsize=(12.2, 10.0))
    heading(fig, "Symbol spelling in selected orthology mappings", "Frozen release: 19 October 2026 | Unique pairwise edges; human-to-source cardinality")
    grid = fig.add_gridspec(2, 2, left=0.09, right=0.94, bottom=0.17, top=0.84,
                           hspace=0.44, wspace=0.30, height_ratios=(1.55, 1))
    for col, graph in enumerate(graphs):
        ax = fig.add_subplot(grid[0, col])
        cards = Counter(r[graph.cardinality_column] for r in graph.rows)
        counts = Counter((r[graph.cardinality_column], graph.spelling(r)) for r in graph.rows)
        wedges, _ = ax.pie([cards[c] for c in CARDINALITIES], radius=0.76,
                           colors=CARD_COLORS, startangle=90, counterclock=False,
                           wedgeprops={"width": 0.27, "edgecolor": "white", "linewidth": 1})
        for wedge, label, card, color in zip(wedges, ("1:1", "1:m", "m:1", "m:m"), CARDINALITIES, CARD_COLORS):
            if cards[card] / graph.edge_count < 0.055:
                continue
            angle = math.radians((wedge.theta1 + wedge.theta2) / 2)
            ax.text(0.625 * math.cos(angle), 0.625 * math.sin(angle), label,
                    ha="center", va="center", fontsize=8, color=INK if color == "#CCBB44" else "white", weight="bold")
        ax.pie([counts[c, s] for c in CARDINALITIES for s in SYMBOLS], radius=1.04,
               colors=list(SYMBOL_COLORS) * len(CARDINALITIES), startangle=90, counterclock=False,
               wedgeprops={"width": 0.255, "edgecolor": "white", "linewidth": 0.7})
        ax.text(0, 0.045, f"{graph.edge_count:,}", ha="center", va="center", fontsize=20, weight="bold")
        ax.text(0, -0.13, "unique edges", ha="center", va="center", fontsize=9, color=GREY)
        ax.set_title(f"{'A' if col == 0 else 'B'}  Human - {graph.title.lower()}", pad=6)
        ax = fig.add_subplot(grid[1, col])
        left = np.zeros(len(CARDINALITIES))
        totals = np.array([cards[c] for c in CARDINALITIES])
        for status, color in zip(SYMBOLS, SYMBOL_COLORS):
            values = np.array([counts[c, status] for c in CARDINALITIES])
            widths = 100 * values / totals
            ax.barh(np.arange(4), widths, left=left, height=0.62, color=color, edgecolor="white", linewidth=0.5)
            for y, width, start, value in zip(range(4), widths, left, values):
                if width >= 15:
                    ax.text(start + width / 2, y, f"{value:,}\n{width:.1f}%", ha="center", va="center", fontsize=8, color="white")
            left += widths
        for y, total in enumerate(totals):
            ax.text(102, y, f"n={total:,}", va="center", fontsize=8)
        ax.set(yticks=range(4), yticklabels=CARD_LABELS, xlim=(0, 122), xticks=(0, 25, 50, 75, 100), xlabel="% of unique edges within cardinality")
        ax.invert_yaxis()
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="x", alpha=0.12)
    fig.legend(handles=legend_patches(CARD_LABELS, CARD_COLORS), loc="upper center", bbox_to_anchor=(0.5, 0.90), ncol=4, frameon=False, fontsize=9, title="Inner ring: final cardinality")
    fig.legend(handles=legend_patches(SYMBOL_LABELS, SYMBOL_COLORS), loc="lower center", bbox_to_anchor=(0.5, 0.093), ncol=3, frameon=False, fontsize=9)
    footer(fig, "Outer ring / bars: compare the recorded symbols after stripping whitespace and converting to uppercase.\nSame spelling is not orthology evidence. Missing either symbol is not classified as different. Counts are edges, not human genes.")
    save(fig, out, "01_symbol_spelling_nested_donuts")


def plot_cardinality(graphs: list[Graph], tables: dict, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 7.1))
    fig.subplots_adjust(left=0.10, right=0.945, bottom=0.25, top=0.79, wspace=0.36)
    heading(fig, "Cardinality: count edges and genes separately", "Selected pairwise graphs only | Each gene is counted once within a class and once in each global total")
    metrics = ("unique_edges", "unique_human_gene_keys", "unique_source_gene_keys")
    colors = ("#68717A", "#228877", "#CC6677")
    labels = ("Unique edges", "Unique human gene keys", "Unique source gene keys")
    for ax, graph, letter in zip(axes, graphs, "AB"):
        rows = {r["final_cardinality"]: r for r in tables["cardinality_counts"] if r["species_pair"] == f"human_{graph.name}"}
        total = rows["All"]
        ax.set_title(f"{letter}  Human - {graph.title.lower()}\n{total['unique_edges']:,} edges | {total['unique_human_gene_keys']:,} human | {total['unique_source_gene_keys']:,} source", fontsize=10.5, loc="left", pad=15)
        for offset, metric, color in zip((-0.24, 0, 0.24), metrics, colors):
            values = [rows[c][metric] for c in CARDINALITIES]
            ax.barh(np.arange(4) + offset, values, height=0.21, color=color)
            for y, value in zip(np.arange(4) + offset, values):
                ax.text(value + 180, y, f"{value:,}", va="center", fontsize=8)
        ax.set(yticks=range(4), yticklabels=CARD_LABELS, xlabel="Count within final human-to-source edge class", xlim=(0, 21000), xticks=(0, 5000, 10000, 15000, 20000))
        ax.invert_yaxis()
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{int(x):,}"))
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="x", alpha=0.12)
    fig.legend(handles=legend_patches(labels, colors), loc="lower center", bbox_to_anchor=(0.5, 0.135), ncol=3, frameon=False, fontsize=9)
    overlap = "; ".join(
        f"{g.title.lower()} {mixed_cardinality_nodes(g, g.human_column):,} human / {mixed_cardinality_nodes(g, g.source_column):,} source"
        for g in graphs
    )
    unresolved = next(r for r in tables["graph_totals"] if r["species_pair"] == "human_fish")["human_keys_without_single_valid_Ensembl_ID"]
    footer(fig, f"Keys in multiple cardinality classes: {overlap}. Class-specific gene counts are not additive.\nGlobal totals are deduplicated, never sums of class counts. Fish counts use canonical fish_key, not exploded Ensembl IDs.\nZebrafish human keys include {unresolved:,} without a single valid Ensembl ID, including a symbol-only unresolved placeholder. No triplet rows are read.")
    save(fig, out, "02_cardinality_edges_and_genes")


def horizontal_counts(ax, values, labels, colors, denominator, title):
    positions = np.arange(len(values))
    ax.barh(positions, values, color=colors, height=0.58)
    for y, value in zip(positions, values):
        percent = 100 * value / denominator
        label = "<0.1%" if 0 < percent < 0.1 else f"{percent:.1f}%"
        ax.text(value + denominator * 0.015, y, f"{value:,}  ({label})", va="center", fontsize=9)
    ax.set(yticks=positions, yticklabels=labels, xlim=(0, denominator * 1.28))
    ax.invert_yaxis()
    ax.tick_params(axis="y", length=0)
    ax.set_title(title, loc="left", fontsize=11, pad=10)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{int(x):,}"))
    ax.grid(axis="x", alpha=0.12)


def plot_fish_evidence(tables: dict, out: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.4))
    fig.subplots_adjust(left=0.17, right=0.95, bottom=0.24, top=0.81, wspace=0.83, hspace=0.85)
    heading(fig, "Zebrafish: selected evidence and feature-ID coverage", "Frozen source-priority policy | Canonical gene keys retained, even when an Ensembl feature cross-reference is unresolved")
    evidence = [r for r in tables["fish_evidence_counts"] if r["final_cardinality"] == "All"]
    n_edges = sum(r["unique_edges"] for r in evidence)
    n_human = sum(r["unique_human_gene_keys"] for r in evidence)
    horizontal_counts(axes[0, 0], [r["unique_edges"] for r in evidence], EVIDENCE_LABELS, EVIDENCE_COLORS, n_edges, f"A  Selection reason | {n_edges:,} unique edges")
    axes[0, 0].set_xlabel("Unique human-fish edges")
    horizontal_counts(axes[0, 1], [r["unique_human_gene_keys"] for r in evidence], EVIDENCE_LABELS, EVIDENCE_COLORS, n_human, f"B  Selection reason | {n_human:,} human keys")
    axes[0, 1].set_xlabel("Unique human gene keys")
    coverage = [r for r in tables["fish_id_coverage"] if r["side"] == "fish" and r["counting_unit"] == "unique_gene_keys"]
    by_status = {r["Ensembl_ID_status"]: r for r in coverage}
    n_fish = coverage[0]["denominator"]
    horizontal_counts(axes[1, 0], [by_status[s]["count"] for s in ID_STATUSES], ID_LABELS, ID_COLORS, n_fish, f"C  Fish Ensembl IDs | {n_fish:,} canonical fish keys")
    axes[1, 0].set_xlabel("Unique canonical fish keys (no ID expansion)")
    usable = [r for r in tables["fish_feature_mapping_counts"] if r["selection_reason"] == "All" and r["coverage_definition"] == "syntax_audited"]
    horizontal_counts(axes[1, 1], [r["unique_edges"] for r in usable], ("Single ID, valid format", "Requires ID review"), (ID_COLORS[0], ID_COLORS[2]), n_edges, f"D  Syntax-audited IDs | {n_edges:,} unique edges")
    axes[1, 1].set_xlabel("Unique human-fish edges")
    cov = {(r["side"], r["counting_unit"], r["Ensembl_ID_status"]): r["count"] for r in tables["fish_id_coverage"]}
    raw_review = next(r["unique_edges"] for r in tables["fish_feature_mapping_counts"] if r["selection_reason"] == "All" and r["coverage_definition"] == "recorded_flags" and not r["usable_for_Ensembl_feature_mapping"])
    reviewed = next(r["unique_edges"] for r in usable if not r["usable_for_Ensembl_feature_mapping"])
    footer(fig, "Curated vs fallback is a selection reason, not an independent vote or a quality ranking. Each human key has one reason.\n"
           f"Fish ID fields: {cov['fish', 'unique_edges', 'missing']:,} edges missing; {cov['fish', 'unique_edges', 'multiple']:,} multiple; {cov['fish', 'unique_edges', 'invalid_syntax']:,} invalid. Human IDs: {cov['human', 'unique_edges', 'missing']:,} edges / {cov['human', 'unique_gene_keys', 'missing']:,} keys missing.\n"
           f"Syntax audit requires review of {reviewed:,} edges vs {raw_review:,} recorded flags: ASS1 contains a ZFIN ID in fish_ensembl_ids; the input is not changed.\n"
           "ID syntax does not confirm current ID existence or biological validity. Canonical fish keys are never expanded into feature IDs.")
    save(fig, out, "03_fish_evidence_and_id_coverage")


def plot_myh6(tables: dict, out: Path) -> None:
    records = tables["myh6_edges"]
    require(len(records) == 2 and {r["species_pair"] for r in records} == {"human_mouse", "human_fish"}, "MYH6 figure expects one direct selected edge per source")
    fig = plt.figure(figsize=(11.8, 5.4))
    heading(fig, "MYH6: two direct pairwise records", "Stable human anchor ENSG00000197616 | Descriptive example from the frozen selected release")
    ax = fig.add_axes((0.045, 0.22, 0.91, 0.60))
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    for row, y, color in zip(records, (0.75, 0.28), ("#68717A", "#4477AA")):
        source = "Mouse" if row["species_pair"] == "human_mouse" else "Zebrafish"
        ax.text(0.02, y + 0.06, f"Human {row['human_symbol']}", weight="bold", fontsize=13, va="center")
        ax.text(0.02, y - 0.03, row["human_ensembl_ids"], fontsize=9, color=GREY, va="center")
        ax.annotate("", xy=(0.64, y), xytext=(0.28, y), arrowprops={"arrowstyle": "-|>", "color": color, "lw": 2})
        ax.text(0.46, y + 0.08, row["final_cardinality"], ha="center", fontsize=11, weight="bold", color=color)
        label = "Preserved frozen mouse pair" if source == "Mouse" else "ZFIN-curated primary"
        ax.text(0.46, y - 0.10, label, ha="center", fontsize=9, color=GREY)
        ax.text(0.68, y + 0.06, f"{source} {row['source_symbol']}", weight="bold", fontsize=13, va="center")
        ids = row["source_key"]
        if source == "Zebrafish":
            ids += "\n" + row["source_ensembl_ids"]
        ax.text(0.68, y - 0.035, ids, fontsize=9, color=GREY, va="top", linespacing=1.6)
    footer(fig, "Fish record: ZFIN_supported=True; archived Ensembl_supported=False; exact pair not listed in Ensembl release 116.\nThese are release-specific annotations, not a universal database failure or evidence of global superiority.\nMYH6 / Myh6 / myh6 have the same spelling after uppercase; that spelling is not orthology evidence. No mouse-fish edge is inferred here.")
    save(fig, out, "04_myh6_pairwise_story")


def checksum(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def export_tables(tables: dict, out: Path) -> None:
    for name, records in tables.items():
        if not records:
            continue
        with (out / f"{name}.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)


def self_test() -> None:
    require([cardinality(a, b) for a, b in ((1, 1), (2, 1), (1, 2), (2, 2))] == list(CARDINALITIES), "Cardinality self-test failed")
    require(symbol_status(" MYH6 ", "myh6") == SYMBOLS[0], "Uppercase comparison failed")
    require(symbol_status("MYH6", "myh7") == SYMBOLS[1], "Different spelling failed")
    require(symbol_status("", "") == SYMBOLS[2], "Empty symbols must not match")
    require(id_status("ENSDARG1;ENSDARG2", "ENSDARG") == "multiple", "Multi-ID coverage failed")
    require(id_status("ENSDARG1;ENSDARG1", "ENSDARG") == "single", "Repeated IDs must be deduplicated")
    require(id_status("", "ENSG") == "missing", "Missing ID coverage failed")
    require(id_status("ZDB-GENE-040912-178", "ENSDARG") == "invalid_syntax", "ZFIN ID must not count as a valid Ensembl ID")
    require(id_status("ENSDARG1;myh6", "ENSDARG") == "invalid_syntax", "Invalid tokens take precedence over multiplicity")
    for fn in (lambda: boolean("unknown"), lambda: cardinality(0, 1)):
        try:
            fn()
        except ValueError:
            continue
        raise AssertionError("Invalid input was not rejected")
    print("Self-tests passed: degrees, spelling, missing symbols, ID cardinality, invalid inputs.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", type=Path, default=DEFAULT_RELEASE)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "docs/figures")
    parser.add_argument("--self-test", action="store_true", help="Run small invariants without reading data or creating figures")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    graphs = [read_graph(args.release_dir / f"human_{name}_final_mapping.tsv", name) for name in ("mouse", "fish")]
    tables = build_tables(graphs)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    style()
    plot_symbols(graphs, args.output_dir)
    plot_cardinality(graphs, tables, args.output_dir)
    plot_fish_evidence(tables, args.output_dir)
    plot_myh6(tables, args.output_dir)
    export_tables(tables, args.output_dir)
    manifest = {
        "description": "Descriptive selected-pairwise orthology mapping figures; no global superiority claim.",
        "credit": "Riley lab / Azad Alizada / CellTalk HHD",
        "frozen_release_date": "2026-10-19",
        "inputs": [{
            "path_relative_to_package": os.path.relpath(g.path.resolve(), ROOT),
            "sha256": checksum(g.path), "size_bytes": g.path.stat().st_size,
            "input_rows": g.raw_rows, "unique_edges": g.edge_count,
        } for g in graphs],
        "script_sha256": checksum(Path(__file__)),
        "runtime": {"python": platform.python_version(), "matplotlib": matplotlib.__version__, "numpy": np.__version__},
        "counting_contract": {
            "mouse_edge": "(human_ensembl_id, mouse_ensembl_id)",
            "fish_edge": "(left_key, fish_key); canonical input keys retained",
            "human_totals": "Unique endpoint keys within each pairwise graph; not a joint-species universe",
            "source_totals": "Unique source endpoint keys; fish feature-ID lists never expanded",
            "cardinality": "Computed from complete selected pairwise graph; stored degrees and labels must agree",
            "symbol_comparison": "Whitespace strip then uppercase; missing either is a separate category",
            "human_symbol_partition": "Missing on any incident edge takes precedence; otherwise all same, all different, or mixed",
            "gene_cardinality_overlap": "Within-class unique gene counts may overlap; All rows deduplicate globally",
            "fish_evidence": "Selection reason, not independent evidence votes; no inference from missing support",
            "feature_usability": "Recorded flags use single-token fields; syntax audit additionally requires ENSG/ENSDARG numeric stable-ID syntax; neither is a biological quality score",
        },
        "validation": {
            "all_rows_selected": "passed", "required_columns_and_canonical_keys": "passed",
            "duplicate_pairs": "exact duplicates deduplicated; conflicting duplicates rejected",
            "stored_degrees_and_cardinality": "passed", "fish_recorded_feature_flags_vs_raw_ID_token_counts": "passed",
            "ID_syntax_issues": len(tables["id_validation_issues"]),
            "one_selection_reason_per_human_fish_key": "passed",
            "consistent_ID_fields_per_canonical_fish_and_human_key": "passed",
            "MYH6_direct_selected_edges": "passed",
        },
        "tables": tables,
    }
    with (args.output_dir / "mapping_audit_counts.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    print(json.dumps({"output_dir": str(args.output_dir.resolve()), "totals": tables["graph_totals"]}, indent=2))


if __name__ == "__main__":
    main()
