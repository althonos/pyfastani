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
__version__ = "0.4.0"

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

# Small addition to the docstring: we want to show a link redirecting to the
# rendered version of the documentation, but this can only work when Python
# is running with docstrings enabled
if __doc__ is not None:
    __doc__ += """See Also:
    An online rendered version of the documentation for this version
    of the library on
    `Read The Docs <https://pytrimal.readthedocs.io/en/v{}/>`_.

    """.format(
        __version__
    )
