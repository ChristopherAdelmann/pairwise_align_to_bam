from collections import defaultdict
import pathlib as pl
import multiprocessing as mp

import pysam as ps

from reference import Reference
from sam_reader import SamReader
from parameter import Parameters

from pairwise_align_to_bam import get_alignment, get_all_alignments
from fasta_reader import FastaReader
from alignment_filter import AlignmentFilter


def process(parameters: Parameters):

    references = FastaReader.parse_references(
        parameters.processing_param.reference_path
    )
    adapter = (
        FastaReader.parse_references(parameters.processing_param.adapter_path)[0]
        if parameters.processing_param.adapter_path
        else None
    )
    header = SamReader.construct_output_header(
        parameters.processing_param.sam_in_path, references
    )

    if adapter:
        parameters.filtering_param.min_offset_from_end = adapter.length

    counts_by_id: defaultdict[str, int] = defaultdict(int)

    total_count = 0
    failed = 0

    with ps.AlignmentFile(
        parameters.processing_param.bam_out_path.as_posix(), "wb", header=header
    ) as align_out, mp.Pool(parameters.processing_param.thread_count) as pool:

        read_generator = SamReader.read_generator(
            parameters.processing_param.sam_in_path,
            references,
            adapter,
            parameters.filtering_param.min_length,
        )

        alignments_chunk: list[ps.AlignedSegment] = []

        if not parameters.filtering_param.all_alignments:
            print("Using optimal alignment")

            for alignment in pool.imap_unordered(get_alignment, read_generator, 1_000):
                total_count += 1

                if total_count % 10_000 == 0:
                    print(
                        f"Total count: {total_count}, Counts by reference: {dict(counts_by_id)}, Failed: {failed}, Fraction failed: {(failed/total_count*100):.2f}%"
                    )
                    for result in alignments_chunk:
                        align_out.write(result)
                    alignments_chunk = []

                if not alignment or not AlignmentFilter.is_valid(
                    alignment, parameters.filtering_param
                ):
                    failed += 1
                    continue

                counts_by_id[alignment.reference_id] += 1
                alignments_chunk.append(alignment.to_sam_record(header))

        else:
            print("Output all alignments filtered by max rel. score difference")

            supplementary_count = 0

            for alignments in pool.imap_unordered(
                get_all_alignments, read_generator, 1_000
            ):
                total_count += 1

                if total_count % 1_000 == 0:
                    print(
                        f"Total count: {total_count}, Failed: {failed}, Fraction failed: {(failed/total_count*100):.2f}%, Supplementary count: {supplementary_count}"
                    )

                    for result in alignments_chunk:
                        align_out.write(result)

                    alignments_chunk = []

                alignments = AlignmentFilter.filtered_alignment_group(
                    alignments, parameters.filtering_param
                )

                if not alignments:
                    failed += 1
                    continue

                supplementary_count += len(alignments) - 1

                alignments_chunk.append(alignments[0].to_sam_record(header))

                for alignment in alignments[1:]:
                    counts_by_id[alignment.reference_id] += 1

                    alignments_chunk.append(alignment.to_sam_record(header, True))

        for result in alignments_chunk:
            align_out.write(result)
