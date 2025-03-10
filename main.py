import sys
import argparse
import pathlib as pl


from process import process
from parameter import *

if __name__ == "__main__":
    import cProfile

    # if check avoids hackery when not profiling
    # Optional; hackery *seems* to work fine even when not profiling, it's just wasteful
    if sys.modules["__main__"].__file__ == cProfile.__file__:
        import pairwise_align_to_bam as pw  # Imports you again (does *not* use cache or execute as __main__)

        globals().update(
            vars(pw)
        )  # Replaces current contents with newly imported stuff
        sys.modules["__main__"] = (
            pw  # Ensures pickle lookups on __main__ find matching version
        )

    parser = argparse.ArgumentParser(
        "PairwiseMap",
        description="Pairwise align ONT reads to a small number of reference sequences",
    )

    parser.add_argument(
        "-i", "--input", type=pl.Path, required=True, help="Input file in SAM format"
    )
    parser.add_argument(
        "-r",
        "--ref",
        type=pl.Path,
        required=True,
        help="Reference file in FASTA format",
    )
    parser.add_argument(
        "-p", "--adpt", type=pl.Path, help="Adapter file in FASTA format"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=pl.Path,
        required=True,
        help="Output file path in BAM format",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=2, help="Number of threads to use"
    )
    parser.add_argument(
        "-l",
        "--minLen",
        type=int,
        default=15,
        help="Minimum length of alignment",
    )
    parser.add_argument(
        "-m",
        "--minMatch",
        type=float,
        default=0.5,
        help="Minimum fraction of bases that need to match",
    )
    parser.add_argument(
        "-a",
        "--all",
        type=bool,
        action=argparse.BooleanOptionalAction,
        help="Output all pairwise alignments and not only the best",
    )
    parser.add_argument(
        "-s",
        "--sup",
        type=float,
        default=1.0,
        help="Max difference from alignments rel. score to optimal alignment rel. score to be output as secondary alignment [0.0-1.0]. Only viable with option --all.",
    )

    args = parser.parse_args()

    processing_params = ProcessingParameters(
        args.threads, args.ref, args.adpt, args.input, args.output
    )

    align_params = AlignmentParameters(5, 8)

    filtering_params = FilteringParameters(
        args.all, args.sup, args.minMatch, args.minLen
    )

    params = Parameters(processing_params, align_params, filtering_params)

    process(params)
