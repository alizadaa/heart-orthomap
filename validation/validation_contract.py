"""API-neutral boundary. No heart-orthomap import is allowed in this module."""

from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple


class AdapterNotReady(RuntimeError):
    pass


@dataclass
class Request:
    X: Any
    source_ids: Tuple[str, ...]
    edges: Tuple[Dict[str, Any], ...]
    method: str = "sum"
    ambiguity: str = "drop"
    adata: Optional[Any] = None
    count_layer: Optional[str] = None
    additional_count_layers: Tuple[str, ...] = ()
    tie_policy: str = "error"
    scores: Optional[Any] = None
    min_score: Optional[float] = None
    feature_key: Optional[str] = None
    feature_type: str = "ensembl"
    output_ids: str = "ensembl"
    strip_ensembl_versions: bool = False


@dataclass
class Result:
    X: Any
    human_ids: Tuple[str, ...]
    adata: Optional[Any] = None
    selected_sources: Optional[Dict[str, Tuple[str, ...]]] = None
    provenance: Optional[Dict[str, Any]] = None
    decisions: Optional[Any] = None
