import pathlib as pl

import biotite.sequence.io.fasta as faio
import biotite.sequence as seq

from reference import Reference


class FastaReader:
    @staticmethod
    def parse_references(path: pl.Path) -> list[Reference]:
        references: list[Reference] = []

        fasta_file = faio.FastaFile.read(path.as_posix())

        for header, sequence in fasta_file.items():
            references.append(Reference(header, seq.NucleotideSequence(sequence)))

        return references
