from dataclasses import dataclass
from functools import cached_property

import biotite.sequence as seq

MAX_MAPQ_SCORE = 60


@dataclass
class Reference:
    reference_id: str
    sequence: seq.NucleotideSequence

    @cached_property
    def length(self) -> int:
        return len(self.sequence)

    @cached_property
    def match_value(self) -> float:
        return MAX_MAPQ_SCORE / self.length
