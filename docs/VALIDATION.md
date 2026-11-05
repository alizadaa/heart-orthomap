# Validation

Local validation date: **5 November 2026**. Code version: **0.1.1**.
These checks validate conversion behavior, not the biological correctness of
every orthology edge or a completed cross-species integration benchmark.

## Executed Checks

| Check | Observed result |
| --- | --- |
| Package regression suite | 24 passed on both tested dependency stacks |
| Portable independent suite | 105 passed on both tested dependency stacks; 79 oracle self-checks passed |
| Independent source execution/oracle suite | 105 passed, 0 failed |
| Independent rebuilt-wheel execution suite | 105 passed, 0 failed; module hashes match source |
| Wheel archive integrity | RECORD hashes verified; no single-cell matrix artifacts bundled |
| Installed wheel, outside the source checkout | Both bundled-map CLI audits and five-strategy synthetic H5AD example passed |
| Bundled mappings | 25,768 mouse edges; 22,361 active fish edges; 645 fish issue rows |
| MYH6 fixture | Curated myh6 link retained despite missing sequence score |
| Sparse stress test | 30,000 cells x 40,000 features, 40,000 nonzeros; exact sparse counts preserved |
| Figure reproducibility | Four PNG/PDF sets; 19 generated artifacts byte-identical across repeat runs |

The independent sparse subprocess used approximately 328 MB peak RSS and
7.78 seconds for conversion on the test machine. Those are measurements of
that sparse fixture, not a portable performance guarantee or production
atlas resource estimate. Direct sparse densification calls were forbidden
during that conversion.

The original test environment was macOS, Python 3.10.19, NumPy 2.2.6, pandas
2.3.3, SciPy 1.15.2 and AnnData 0.11.4. A second isolated environment tested
Python 3.12.14, NumPy 2.3.5, pandas 3.1.1, SciPy 1.18.1 and AnnData 0.13.4.
The 0.1.1 correction passed the package and portable independent suites in
both environments. Original wheel/stress measurements above describe the
0.1.0 review; they are not newly measured 0.1.1 wheel performance claims.

The first Linux/Python 3.10-3.12 GitHub CI run exposed an AnnData compatibility
failure. Version 0.1.1 fixes that failure. The corrected
[Linux matrix passed on all three Python versions](https://github.com/alizadaa/heart-orthomap/actions/runs/42110829202).
Local editable-install checks remain distinct from clean Linux installations.

## Review Findings Fixed

1. Source-ID collision quarantine formerly removed rows before checking target
   collisions, potentially hiding overlapping ambiguity. Both collision sets
   are now computed on the same original eligible graph and removed together.
2. External numeric-string expression scores were validated as numeric but
   merged as strings. They are now copied and normalized before thresholding.
3. AnnData 0.13 includes `.X` as `None` in layer iteration. The conversion's
   named-layer audit now excludes that sentinel, without changing counts.
   New regression tests cover explicit raw `.X`, retained count layers and
   H5AD write/read.

The portable fixture also respects current AnnData's CSR/CSC-only storage:
COO inputs are checked for explicit rejection when the AnnData constructor
does not support them. Older AnnData still exercises COO conversion. The
zero-library fixture builds a temporary LIL matrix locally and assigns only
its CSR result. The independent arithmetic oracle and golden values did not
change; unsupported storage is not silently coerced by the adapter.

Both cases have permanent package regression tests. The independent reviewer
did not edit the implementation; fixes were rerun against source and the built
wheel. All five strategies also passed synthetic H5AD write/read roundtrips,
including observation order, requested layers and cell-aligned embeddings.

## Reproduce

From an installed checkout:

```bash
python -m unittest discover -s tests -v
python -B validation/run_validation.py --self-check
python -B validation/run_validation.py --installed
python examples/convert_toy.py --out toy-results
heart-orthomap audit --species mouse --out mouse-audit
heart-orthomap audit --species zebrafish --out fish-audit
```

Use new output directories. Original full mapping tables are needed only to
rebuild the compact frozen data or regenerate figures from their upstream
graph; ordinary conversion and audit commands use the bundled frozen data.

The [portable independent harness](../validation/README.md) includes its own
integer-accumulation oracle, synthetic fixtures and both bug reproducers. It
tests the installed public API rather than injecting a source path or repairing
package outputs. Generated reports and synthetic H5ADs are ignored by Git.

## Limits

No patient counts or real cell barcodes are included in validation fixtures.
Count conservation does not establish the best biological paralog-selection
policy. Native QC/library-size provenance, annotations, independent benchmark
labels, upstream data terms, other operating systems and future dependency
versions require their own checks. Historical integration matrices are not
retroactively certified by this package's tests.
