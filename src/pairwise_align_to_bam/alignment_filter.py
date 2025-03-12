from pairwise_align_to_bam.read import AlignedRead
from pairwise_align_to_bam.parameter import FilteringParameters


class AlignmentFilter:

    @staticmethod
    def is_valid(alignment: AlignedRead, parameters: FilteringParameters) -> bool:
        ref_start_pos: int = alignment.alignment.trace[0][0]
        offset_from_ref_end = len(alignment.alignment.sequences[0]) - ref_start_pos

        return (
            alignment.rel_align_score >= parameters.min_rel_score
            and alignment.aligned_len() >= parameters.min_length
            and offset_from_ref_end >= parameters.min_offset_from_end
        )

    @staticmethod
    def filtered_alignment_group(
        alignments: list[AlignedRead], parameters: FilteringParameters
    ) -> list[AlignedRead]:
        if not alignments:
            return []

        alignments.sort(reverse=True)

        max_rel_align_score = alignments[0].rel_align_score

        valid_alignments: list[AlignedRead] = []

        for alignment in alignments:
            if not AlignmentFilter.is_valid(alignment, parameters):
                continue

            if (
                not (max_rel_align_score - alignment.rel_align_score)
                <= parameters.max_diff_from_optimal
            ):
                break

            valid_alignments.append(alignment)

        return valid_alignments
