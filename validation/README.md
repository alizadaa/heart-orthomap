# Independent heart-orthomap validation

Portable reviewer tests for the public heart-orthomap 0.1.0 API. The suite uses
the normally installed `heart_orthomap`; it does not search a checkout, modify
import paths, patch the package, or query the network. It can live in
`tests/reviewer/` or any other writable directory.

## Run

Use Python 3.10+ in an environment with heart-orthomap and its dependencies
(NumPy, SciPy, pandas, AnnData/h5py). Install the wheel to be reviewed into that
environment first, for example:

```sh
python -m pip install ./heart_orthomap-0.1.0-py3-none-any.whl
```

From this directory:

```sh
python -B run_validation.py --self-check
python -B run_validation.py --installed
python -B reproduce_findings.py numeric-string-scores
python -B reproduce_findings.py overlapping-collisions
```

Expected on the reviewed 0.1.0 API: **79 oracle self-checks** and **105 installed
package checks**, all passing; both regression commands exit zero. A missing
package/dependency or a failed assertion is an error, not a skipped success.
Run without `-O` or `PYTHONOPTIMIZE`; optimized runs are rejected explicitly.
The self-check validates the oracle/harness, not the installed package.

Optional separate-process sparse scale check on macOS/Linux:

```sh
python -B run_sparse_stress.py
```

It maps a synthetic 30,000 x 40,000 matrix with 40,000 nonzeros, forbids direct
sparse densification during conversion, and checks exact output and peak RSS.
The process must stay below a deliberately generous 2 GiB RSS ceiling; measured
RSS depends on platform/import overhead. No dense 9.6 GB matrix is allocated.

## Inputs and outputs

- Distributed inputs are three small JSON files containing entirely synthetic
  IDs/counts and hand-calculated expectations, plus Python fixture generators.
- No real count data, cell barcodes, H5AD files, curated mapping copies, reports,
  binaries, environment files, or original workspace paths are distributed.
- The bundled-species sanity check reads mapping resources from the installed
  package only. It does not include or download copies of those resources.
- Runs create `reports/` beside these scripts, including JSON results, generated
  synthetic tabular fixtures and five synthetic H5AD roundtrips. This directory
  is ignored by Git and must not be included when redistributing the harness.
- Reports record runtime environment, imported-package code hashes before/after,
  fixture hashes, individual outcomes and exceptions. Nonzero failures or code
  changes produce a nonzero exit. A report is not a universal release approval.

## Coverage and independence

The golden count fixtures and adversarial graph design were made before
implementation inspection. The adapter translates names/return structures only;
it never repairs count arithmetic, deduplication, ambiguity or score selection.
The reference oracle uses Python integer accumulation rather than the package's
projection algorithm. `METHODS.md` states the executed policies and limitations.

Checks include all four human-perspective types; repeated Cartesian rows;
ambiguity drop/error; all five conversion strategies; ties/missing scores;
stable IDs and symbol collisions; invalid/unsafe numbers; sparse formats;
permutations; immutability; metadata/layer alignment; zero mapping; matched-group
expression scoring; native readers; both fixed bug reproducers; and H5AD I/O.
Bundled-data sanity assertions are specific to the reviewed release and should
be explicitly revised for a new curated release, never silently removed.
