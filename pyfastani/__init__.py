from . import _fastani
from ._fastani import Sketch, Mapper, Hit, MAX_KMER_SIZE

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "MIT"
__version__ = "0.2.0"

__all__ = ["Sketch", "Mapper", "Hit", "MAX_KMER_SIZE"]
__doc__ = _fastani.__doc__
