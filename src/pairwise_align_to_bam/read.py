from dataclasses import dataclass
from functools import cached_property, total_ordering
import math
from venv import logger

import biotite.sequence as seq
import biotite.sequence.align as align
import pysam as ps
import numpy as np


@dataclass
class Read:
    read_id: str
    sequence: seq.NucleotideSequence
    quality_string: str
    tags: str

    @cached_property
    def length(self):
        return len(self.sequence)

    def relative_alignment_score(self, alignment_score: int, match_score: int) -> float:
        return alignment_score / (self.length * match_score)


@dataclass
@total_ordering
class AlignedRead(Read):
    reference_id: str
    alignment: align.Alignment
    rel_align_score: float

    @classmethod
    def from_read(
        cls,
        read: Read,
        reference_id: str,
        alignment: align.Alignment,
        rel_align_score: float,
    ):
        return cls(
            read_id=read.read_id,
            sequence=read.sequence,
            quality_string=read.quality_string,
            tags=read.tags,
            reference_id=reference_id,
            alignment=alignment,
            rel_align_score=rel_align_score,
        )

    def aligned_len(self) -> int:
        return len(self.alignment.trace)

    def to_sam_record(
        self, sam_header: ps.AlignmentHeader, is_supplementary=False
    ) -> ps.AlignedSegment:
        alignment_start_position = self.alignment.trace[0][0] + 1

        cigar_string: str | np.ndarray = align.write_alignment_to_cigar(
            self.alignment,
            distinguish_matches=True,
            include_terminal_gaps=True,
            as_string=True,
        )

        assert type(cigar_string) == str

        flag = "2048" if is_supplementary else "0"

        # Probability values for modified bases need to be handled separately because pysam cannot currently
        # convert it from string to alignment segment (https://github.com/pysam-developers/pysam/issues/1123)
        ml_values: list[int] | None = None
        ml_tag = "ML:B:C"

        if ml_tag in self.tags:
            split_tags = self.tags.split("\t")

            cleaned_tags: list[str] = []

            for split_tag in split_tags:
                if ml_tag in split_tag:
                    ml_values = [
                        int(val) for val in split_tag.split(":")[-1].split(",")[1:]
                    ]
                else:
                    cleaned_tags.append(split_tag)

            self.tags = "\t".join(cleaned_tags)

        sam_string = "\t".join(
            [
                self.read_id,
                flag,
                self.reference_id,
                str(alignment_start_position),
                "0",
                cigar_string,
                "*",
                "0",
                "0",
                str(self.sequence),
                self.quality_string,
                f"sc:i:{self.alignment.score}",
                f"rs:f:{self.rel_align_score}",
                self.tags,
            ]
        )

        aligned_segment = ps.AlignedSegment.fromstring(sam_string, sam_header)

        if ml_values and len(ml_values) > 0:
            aligned_segment.set_tag("ML", ml_values)

        return aligned_segment

    def __eq__(self, other: "AlignedRead") -> bool:
        return math.isclose(self.rel_align_score, other.rel_align_score)

    def __lt__(self, other: "AlignedRead") -> bool:
        return self.rel_align_score < other.rel_align_score
