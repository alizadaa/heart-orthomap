#!/usr/bin/env python3
"""Run separately; direct sparse densification is forbidden during conversion."""

import json
from pathlib import Path
import sys
import time

sys.dont_write_bytecode = True
if not __debug__:
    raise RuntimeError("Assertions are required; run without -O or PYTHONOPTIMIZE")
try:
    import resource
except ImportError:
    raise SystemExit("The optional RSS stress check requires a Unix-like Python")

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

import package_adapter

ROOT = Path(__file__).resolve().parent


def main():
    api = package_adapter.api()
    hashes_before = package_adapter.source_hashes()
    rows, cols = 30000, 40000
    source_ids = ["STRESS_S%05d" % i for i in range(cols)]
    human_ids = ["STRESS_H%05d" % i for i in range(cols)]
    values = np.arange(cols, dtype=np.int64) % 97 + 1
    matrix = sparse.coo_matrix((values, (np.arange(cols) % rows, np.arange(cols))), shape=(rows, cols)).tocsr()
    data = ad.AnnData(matrix, var=pd.DataFrame(index=source_ids))
    mapping = api.MappingTable(pd.DataFrame({
        "source_id": source_ids, "target_id": human_ids,
        "source_symbol": source_ids, "target_symbol": human_ids,
    }), "mouse")
    baseline = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    originals = []

    def no_dense(*args, **kwargs):
        raise AssertionError("Sparse matrix densification attempted during mapping")

    for cls in (sparse.csr_matrix, sparse.csc_matrix, sparse.coo_matrix):
        for name in ("toarray", "todense"):
            originals.append((cls, name, getattr(cls, name)))
            setattr(cls, name, no_dense)
    start = time.monotonic()
    try:
        result = api.convert_anndata(data, mapping, strategy="sum", counts_layer=None)
    finally:
        for cls, name, original in originals:
            setattr(cls, name, original)
    elapsed = time.monotonic() - start
    output = result.adata.X
    assert sparse.issparse(output) and output.shape == (rows, cols) and output.nnz == cols
    actual = output.copy()
    actual.sort_indices()
    np.testing.assert_array_equal(actual.data, matrix.data)
    np.testing.assert_array_equal(actual.indices, matrix.indices)
    np.testing.assert_array_equal(actual.indptr, matrix.indptr)
    assert tuple(result.adata.var_names) == tuple(human_ids)
    assert int(output.sum()) == int(values.sum())
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    units = 1 if sys.platform == "darwin" else 1024
    hashes_after = package_adapter.source_hashes()
    report = {"status": "passed", "shape": [rows, cols], "nnz": cols,
              "dense_float64_bytes_not_allocated": rows * cols * 8,
              "direct_sparse_densification_forbidden": True,
              "baseline_peak_rss_bytes": baseline * units, "peak_rss_bytes": peak * units,
              "elapsed_seconds": elapsed, "source_hashes_before": hashes_before,
              "source_hashes_after": hashes_after, "source_unchanged": hashes_before == hashes_after,
              "note": "Separate subprocess; observed RSS, not a portable hard memory cap."}
    assert report["peak_rss_bytes"] < 2 * 1024**3, "Unexpected >2 GiB process RSS"
    (ROOT / "reports").mkdir(exist_ok=True)
    (ROOT / "reports" / "sparse_stress.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
