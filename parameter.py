from dataclasses import dataclass
import pathlib as pl


@dataclass
class ProcessingParameters:
    thread_count: int
    reference_path: pl.Path
    adapter_path: pl.Path | None
    sam_in_path: pl.Path
    bam_out_path: pl.Path


@dataclass
class AlignmentParameters:
    gap_open_penalty: int
    gap_extension_penalty: int


@dataclass
class FilteringParameters:
    all_alignments: bool
    max_diff_from_optimal: float
    min_rel_score: float
    min_length: int
    min_offset_from_end: int = 0


@dataclass
class Parameters:
    processing_param: ProcessingParameters
    alignment_param: AlignmentParameters
    filtering_param: FilteringParameters
