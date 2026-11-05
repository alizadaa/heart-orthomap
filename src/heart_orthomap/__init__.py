"""Heart-focused orthology conversion with explicit count and ambiguity policies."""
from .mapping import MappingTable, load_mapping, load_bundled_mapping
from .convert import ConversionResult, convert_anndata
from .expression import expression_similarity_scores

__version__ = "0.1.0"
__all__ = ["MappingTable", "load_mapping", "load_bundled_mapping", "ConversionResult",
           "convert_anndata", "expression_similarity_scores"]
