import biotite.sequence.align as align
import numpy as np

from reference import Reference
from read import Read, AlignedRead


def get_effective_length(trace: np.ndarray) -> int:
    start = 0
    while start < len(trace) and trace[start, 1] == -1:
        start += 1

    # Find end index (last row from the end with second element not -1)
    end = len(trace) - 1
    while end >= 0 and trace[end, 1] == -1:
        end -= 1

    # Slice the array to keep only the rows between start and end (inclusive)
    return np.count_nonzero(trace[start : end + 1][:, 1] != -1)


def get_relative_alignment_score(alignment: align.Alignment, match_score: int) -> float:
    alignment_score = alignment.score or 0

    effective_len = len(alignment.sequences[1])  # get_effective_length(alignment.trace)

    return alignment_score / (effective_len * match_score)


def align_read(read: Read, reference: Reference) -> AlignedRead | None:
    matrix = align.SubstitutionMatrix.std_nucleotide_matrix()
    match_score = matrix.score_matrix()[0][0]

    alignment: align.Alignment = align.align_optimal(
        reference.sequence,
        read.sequence,
        matrix,
        gap_penalty=(-5, -8),
        terminal_penalty=True,
        local=True,
        max_number=1,
    )[0]

    if not alignment.score:
        return None

    relative_alignment_score = get_relative_alignment_score(alignment, match_score)

    return AlignedRead.from_read(
        read,
        reference.reference_id,
        alignment,
        relative_alignment_score,
    )


def is_adapter(adapter: Reference, read: Read, min_identity: float) -> bool:

    alignment: AlignedRead | None = align_read(read, adapter)

    if not alignment:
        return False

    max_distance_from_start = ((1 - min_identity) * adapter.length) * 4

    align_seq_start = alignment.alignment.trace[0][1]

    return (
        align.get_sequence_identity(alignment.alignment, "all") >= min_identity
        or max_distance_from_start > align_seq_start
    )


def get_alignment(
    input: tuple[list[Reference], Reference | None, Read]
) -> AlignedRead | None:
    references, adapter, read = input

    # TODO: Make this a command line parameter
    if adapter and is_adapter(adapter, read, 0.95):
        return None

    best_alignment: None | AlignedRead = None

    for reference in references:
        alignment: AlignedRead | None = align_read(read, reference)

        if not alignment:
            continue

        if not best_alignment or alignment > best_alignment:
            best_alignment = alignment

    return best_alignment


def get_all_alignments(
    input: tuple[list[Reference], Reference | None, Read]
) -> list[AlignedRead]:
    references, adapter, read = input

    # TODO: Make this a command line parameter
    if adapter and is_adapter(adapter, read, 0.95):
        return []

    alignments: list[AlignedRead] = []

    for reference in references:
        alignment = align_read(read, reference)

        if alignment:
            alignments.append(alignment)

    return alignments
