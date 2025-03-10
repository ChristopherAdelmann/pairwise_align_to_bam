from typing import Generator
import pathlib as pl

import pysam as ps
import biotite.sequence as seq

from reference import Reference
from read import Read


class SamReader:

    @staticmethod
    def construct_output_header(
        input_sam_path: pl.Path, references: list[Reference]
    ) -> ps.AlignmentHeader:
        header_string = ""

        with open(input_sam_path, "r") as file:
            for line in file:
                if not line.startswith("@"):
                    break

                header_string += line

        for reference in references:
            header_string += "\t".join(
                ["@SQ", f"SN:{reference.reference_id}", f"LN:{len(reference.sequence)}"]
            )
            header_string += "\n"

        return ps.AlignmentHeader.from_text(header_string)

    @staticmethod
    def read_generator(
        path: pl.Path,
        references: list[Reference],
        adapter: Reference | None,
        min_length: int,
    ) -> Generator[tuple[list[Reference], Reference | None, Read], None, None]:
        with open(path, "r", buffering=8192) as file:
            for line in file:
                if line.startswith("@"):
                    continue

                tokens = line.split("\t", 11)

                if len(tokens[9]) < min_length:
                    continue

                yield (
                    references,
                    adapter,
                    Read(
                        tokens[0],
                        seq.NucleotideSequence(tokens[9]),
                        tokens[10],
                        tokens[11] if len(tokens) > 11 else "",
                    ),
                )
