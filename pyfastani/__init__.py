from . import _fastani
from ._fastani import (
    Sketch,
    Mapper,
    Hit,
    Minimizers,
    MinimizerInfo,
    MinimizerIndex,
    Position,
    MAX_KMER_SIZE
)

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "MIT"
__version__ = "0.3.0"

__doc__ = _fastani.__doc__
__all__ = [
    "Sketch",
    "Mapper",
    "Hit",
    "Minimizers",
    "MinimizerInfo",
    "MinimizerIndex",
    "Position",
    "MAX_KMER_SIZE"
]
