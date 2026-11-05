# Executed validation contract

## Graph and count semantics

Edges are distinct human/source stable-ID pairs, not third-species triplet
rows. The four disconnected motifs represent human-perspective One2One,
One2Many (one human, several sources), Many2One (several humans, one source),
and Many2Many. Synthetic triplet repetition must not multiply count weights.

A source linked to several humans is explicitly dropped or rejected; no method
copies its full count to several targets. The default strict One2One strategy
is tested separately from requested non-strict strategies. Sum combines all
eligible distinct source columns per human. Representative strategies choose
one global source per human, never a separate per-cell maximum.

`best_sequence` uses the arithmetic mean of both directional percentage
identities, with no current-to-archived fallback. `best_expression` uses explicit
matched-group expression correlations, not abundance. The separate legacy
`highest_total_expression` sensitivity strategy uses total raw counts across
all supplied observations. Tied or missing comparison scores trigger the
explicit error/drop policy; there is no arbitrary stable-ID winner. Singleton
candidates need no comparison unless a threshold is requested.

Stable IDs govern identity; shared symbols never merge distinct genes.
Native-reader fixtures use artificial, namespace-valid 99-prefixed Ensembl
IDs. Canonical in-memory fixtures use short synthetic IDs. Conflicting pair
metadata is rejected. Source/target canonical-ID collisions are detected on
the same eligible baseline and quarantined atomically.

## Arithmetic and independent oracle

Stored raw counts must be nonnegative finite integers, with entries strictly
below 2^53 under this API. Whole-number floats within the safe range are
accepted; fractional, nonfinite, negative, bool/complex and unsafe values are
rejected explicitly. int32 input promotes safely. Sums can exceed 2^53 when
exactly represented as int64; unsafe totals must not wrap or round. The oracle
accumulates Python integers rather than sharing the package's sparse matrix
projection. Hand-calculated golden matrices check destinations, not just totals.

Inputs must have unique nonempty feature identities. Explicit version stripping
must reject resulting collisions. Canonical CSR/CSC/COO counts are tested;
noncanonical sparse duplicate entries are rejected by the confirmed API.
No surviving target is a clear ValueError, not an invalid matrix; zero
observations and mapped all-zero genes are supported.

## Alignment, scoring and I/O

Original cell order, observation columns, cell embeddings and graphs are
preserved. An audit library-size column may be added. Output X/counts are raw
mapped counts. Requested extra count layers use the same mapping/selection
plan, even when their own abundance rankings disagree. Unrequested layers and
source raw/varm/varp/uns are dropped and reported. Inputs remain unchanged.
All five strategies undergo synthetic H5AD write/read roundtrips.

Expression similarity is independently calculated with scalar sums: normalize
each full-native cell library to CP10K, average within matched supported groups,
apply log1p, then compute Pearson r across equally weighted groups. Group-Series
indices must exactly match the ordered observation axis. Constant profiles
produce NaN; unsupported groups, zero libraries and invalid labels are rejected.

## Interpretation

Passing verifies the exercised arithmetic/API contract, not biological validity
of every orthology edge or superiority to another method. Highest-total
selection remains a legacy sensitivity choice, not formal biological evidence.
No patient or experimental single-cell data are distributed. Resource and
runtime measurements are environment-specific; arbitrary type-filter controls
are not exposed by this API, and cross-version/platform coverage is not claimed.
