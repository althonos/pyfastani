import typing
from typing import Dict, Generic, List, Union, Iterable, Tuple, TypeVar

N = TypeVar("N")
Seq = Union[str, bytes, bytearray, memoryview]

MAX_KMER_SIZE: int

class _Parameterized:
    @property
    def k(self) -> int: ...
    @property
    def window_size(self) -> int: ...
    @property
    def fragment_length(self) -> int: ...
    @property
    def minimum_fraction(self) -> float: ...
    @property
    def percentage_identity(self) -> float: ...
    @property
    def p_value(self) -> float: ...
    @property
    def protein(self) -> bool: ...

class Sketch(Generic[N], _Parameterized):
    def __init__(
        self,
        *,
        k: int = 16,
        fragment_length: int = 3000,
        minimum_fraction: float = 0.2,
        p_value: float = 1e-03,
        percentage_identity: float = 80.0,
        reference_size: int = 5_000_000,
        protein: bool = False,
    ) -> None: ...
    def __getstate__(self) -> Dict[str, object]: ...
    def __setstate__(self, state: Dict[str, object]) -> None: ...
    @property
    def occurences_threshold(self) -> int: ...
    @property
    def names(self) -> List[N]: ...
    @property
    def minimizers(self) -> List[MinimizerInfo]: ...
    def add_draft(self, name: N, contigs: Iterable[Seq]) -> Sketch[N]: ...
    def add_genome(self, name: N, sequence: Seq) -> Sketch[N]: ...
    def clear(self) -> Sketch[N]: ...


class Mapper(Generic[N], _Parameterized):
    def __getstate__(self) -> Dict[str, object]: ...
    def __setstate__(self, state: Dict[str, object]) -> None: ...
    @property
    def minimizers(self) -> List[MinimizerInfo]: ...
    def query_draft(self, contigs: Iterable[Seq]) -> List[Hit[N]]: ...
    def query_genome(self, sequence: Seq) -> List[Hit[N]]: ...


class Hit(Generic[N]):
    def __init__(self, name: N, identity: float, matches: int, fragments: int): ...
    def __repr__(self) -> str: ...
    def __eq__(self, other: object) -> bool: ...
    def __getstate__(self) -> Tuple[N, float, int, int]: ...
    def __setstate__(self, state: Tuple[N, float, int, int]) -> None: ...
    @property
    def name(self) -> N: ...
    @property
    def matches(self) -> int: ...
    @property
    def fragments(self) -> int: ...
    @property
    def identity(self) -> float: ...


class MinimizerInfo:
    def __init__(self, hash: int, sequence_id: int, window_position: int): ...
    def __repr__(self) -> str: ...
    def __eq__(self, other: object) -> bool: ...
    def __getstate__(self) -> Tuple[int, int, int]: ...
    def __setstate__(self, state: Tuple[int, int, int]) -> None: ...
    @property
    def hash(self) -> int: ...
    @property
    def sequence_id(self) -> int: ...
    @property
    def window_position(self) -> int: ...
