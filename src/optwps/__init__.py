from .optwps import WPS, main
from .utils import exopen, is_soft_clipped, ref_aln_length

import importlib.metadata

try:
    __version__ = importlib.metadata.version(__name__)
except importlib.metadata.PackageNotFoundError:
    __version__ = "0.0.0"  # Fallback for development mode

__all__ = ["WPS", "exopen", "is_soft_clipped", "main", "ref_aln_length"]
