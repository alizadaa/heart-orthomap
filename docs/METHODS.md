# Methods

## Scope and Frozen Data

Version 0.1 converts mouse or zebrafish RNA gene features to human using a
frozen, human-anchored curated graph. It does not perform cross-species
integration, infer cell lineage, change cell annotations, or project ATAC peaks.
Chicken and rat are not part of this first packaged data release.

The human-mouse graph preserves the original Ensembl stable-ID edge set.
For each human gene, direct ZFIN human-zebrafish orthologs take priority; only
when that human gene has no ZFIN link are its archived Ensembl links retained.
Alliance 9.2, HCOP and NCBI provide supporting annotations in the source release,
not independent votes or automatic replacement edges. Ensembl release 116 adds
exact-pair scores without changing the selected graph. Sequence identity is
not a calibrated orthology probability.

`MYH6` to `myh6` illustrates why curated zebrafish evidence matters: this pair
is selected by ZFIN even though it has no exact-pair identity score in the
frozen current-Ensembl export. A missing score is not a zero score and does
not veto an unambiguous curated pair. This example does not prove that every
curated assignment is correct or that the resource is universally superior.

## Pairwise Graphs, Not Repeated Triplets

The final selected three-species TSV contains 25,590 triplet rows, but only
18,202 distinct human-mouse and 20,891 distinct human-fish pairs. Its human-mouse
coverage is narrower than the complete human-mouse graph (25,768 pairs); the
complete human-fish graph has 23,006 pairs. Therefore the bundled data are built
from the complete pairwise final tables, not the selected triplet view.

The annotated triplet reader honors `human_to_mouse_final_mapping_selected`
or `human_to_fish_final_mapping_selected` independently. It does not use the
fish selection flag to remove otherwise valid human-mouse pairs. Repeated
stable-ID pairs are collapsed once; conflicting pair metadata raise an error.
The stored complete-graph degrees are retained rather than recomputed from a
three-species subset or an input dataset's currently measured genes.

## Direction and Cardinality

For an edge between human gene h and source gene s, let dH(h) be the number
of source neighbors and dS(s) the number of human neighbors in the selected
complete graph:

| Human-perspective class | dH | dS |
| --- | --- | --- |
| One2One | 1 | 1 |
| One2Many | >1 | 1 |
| Many2One | 1 | >1 |
| Many2Many | >1 | >1 |

These are edge-local classes, not labels for entire connected components.
In mouse/fish-to-human conversion, human One2Many is the situation in which
several source genes can collapse into one human feature. A source gene with
several human neighbors is a different problem: expanding it would duplicate
its count vector. This release rejects such sources or drops them under an
explicit policy, consistently across all strategies. It does not claim a
globally optimal matching, fractional allocation, or automatic resolution of
Many2Many relationships.

## Identifier Safety

Stable Ensembl IDs are preferred. Exact source symbols may be used only when
they identify one stable source gene; no uppercasing or alias guessing occurs.
Optional Ensembl version stripping is explicit, and resulting collisions fail.
Missing or multiple feature IDs are quarantined, not expanded into extra genes.

The source fish graph flags 638 selected rows with unresolved/multiple
feature IDs. One further row (`ASS1-ass1`) contains a ZFIN ID in its Ensembl
field, which the package quarantines as an invalid identifier namespace.
An additional six rows share three Ensembl IDs between different
canonical genes (`cenpq/mrpl39`, `pdcd7/ubap1lb`, `fasn/ccdc57`). These six rows
are also quarantined. The packaged resolved fish projection consequently has
22,361 edges; this is distinct from the 23,006-edge curated canonical graph.
The original tables are not edited. The unresolved-ID table accompanies the
bundle for review and future versioned corrections.

## Count Projection

Let X be a cells-by-native-features raw-count matrix, and P a binary
native-features-by-human-features projection. The output is Y = X P.
Every row of P has at most one nonzero entry, so one native feature is never
counted twice across human targets. In Sum mode, a column of P may have
multiple nonzero entries. Representative modes allow at most one.

The tool validates finite, nonnegative integer counts, checks arithmetic
bounds, and uses int64 sparse multiplication. It does not sum normalized data,
round fractional data to counts, average p-values, or normalize the output.
For every cell, the output count sum must exactly equal the input count sum
over the selected source columns. It does not promise equality to the full
input library, because unmapped/ambiguous/unselected genes are intentionally
excluded. Source count arrays are not modified.

## Strategies

1. **One2One:** retain only final reciprocal One2One edges. This is the
   conservative default and a useful baseline.
2. **Sum:** retain all eligible unambiguous source features and sum their raw
   per-cell counts into each human target.
3. **Best sequence:** among source candidates present in the input, select
   the unique highest arithmetic mean of the two exact-pair percentage
   identities. The default source is Ensembl 116. An archived-score mode can
   be explicitly requested where that source is available. There is no silent
   current-to-archived fallback, fuzzy column-name detection, or invented
   score for a curated pair.
4. **Best expression:** select the unique highest explicitly supplied
   cross-species expression correlation. The helper computes Pearson r over
   `log1p(mean native CP10K)` profiles in at least three adequately supported,
   explicitly matched biological groups. It normalizes each cell using all
   supplied native features before group averaging. Cells are pooled within
   each group and groups receive equal weight; donors are not automatically
   balanced. Constant profiles yield missing scores. Unpaired individual
   cells are never correlated across species. Use separate training/reference
   cells, not benchmark holdout labels.
5. **Highest total expression:** a separately named legacy-style sensitivity,
   choosing the largest raw-count total over the supplied input cells. This
   is abundance selection, not expression similarity. It is dataset/sample
   dependent; pass the intended population explicitly. Selection is once per
   human target, not switched per cell.

When only one candidate is present there is no ranking choice, so it is kept
without requiring a score unless a minimum score was explicitly requested.
With several candidates, missing comparison scores or tied maxima raise an
error by default. The optional drop policy excludes that target and records
why. Minimum-score thresholds apply only to representative-selection modes.
Every selection/drop appears in the decisions table.

## Historical Compatibility and Fixes

The older expression mode used sample-wide total counts, not similarity to a
human expression profile. The older Python sequence reader could match the
same identity column twice; the new reader uses both exact column names.
An older R Sum path could repeat symbol pairs after a metadata join; the new
projection deduplicates stable-ID pairs. Older scripts could copy one source
vector into several human features; this release never does so. The project-
specific RHD-to-RHCE alias is not imported as biological orthology.

Consequently, the new implementation is intentionally not a byte-for-byte
replication of every historical matrix. Regression fixtures establish count
semantics, not which script version generated every historical integration.

## Preserved and Rebuilt Objects

The returned AnnData has identical cell IDs/order and original `.obs` columns,
plus an input-library-size accounting column. Cell-aligned `.obsm` and `.obsp`
are preserved, with provenance that they are native embeddings/graphs, not
recomputed on the human feature matrix. Feature metadata are rebuilt with human
IDs, symbols, selected source IDs and source-feature multiplicities. Gene-axis
objects and stale feature analyses are not silently relabeled. Only explicitly
requested additional raw count layers are projected through the same P.

Keep the original native object. For downstream comparisons, use original
full-native library sizes rather than re-normalizing over different mapped or
HVG-restricted universes. The tool cannot establish that an input already has
its full native gene axis; the expression-scoring helper requires that scope
as an explicit input assumption. All downstream DE and biological validation
remain separate from orthology conversion.
