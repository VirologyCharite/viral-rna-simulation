import argparse

from viral_rna_simulation.simulate import run


def parse_args() -> argparse.Namespace:
    """
    Make an argument parser and use it to parse the command line.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Simulate the polymerization of viral RNA within cells and the subsequent "
            "RNA sample preparation, sequencing, and alignment to the viral genome."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--cells",
        type=int,
        default=1,
        metavar="N",
        help=(
            "The number of cells to simulate. This is purely for speed-up so multiple "
            "cells can be run in parallel. All cells are seeded with one copy of the "
            "same (+) RNA molecule."
        ),
    )

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument(
        "--genome-length",
        type=int,
        metavar="N",
        help=(
            "Specify the viral RNA genome length. A random genome of the given length "
            "will be used. Incompatible with --genome."
        ),
    )

    group.add_argument(
        "--genome",
        metavar="ACGT...",
        help="Specify the infecting viral genome. Incompatible with --genome.",
    )

    parser.add_argument(
        "--steps",
        type=int,
        default=1000,
        metavar="N",
        help=(
            "The number of replication steps to simulate. In each replication step, a "
            "random RNA from each cell will be chosen to replicate. The replication "
            "will create the reverse complement sequence with the (+/-) sense flipped. "
            "The new molecule will be added to the RNA population within the cell."
        ),
    )

    parser.add_argument(
        "--mutation-rate",
        type=float,
        default=0.001,
        metavar="N",
        help=(
            "The per-nucleotide mutation (polymerase misincorporation) rate, applied "
            "during RNA molecule replication."
        ),
    )

    parser.add_argument(
        "--ratio",
        type=int,
        default=1,
        metavar="N",
        help=(
            "The number of (+) RNA molecules to make from each (-) molecule. In "
            "coronaviruses this ratio is often in the range of 10 to 100."
        ),
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run(
        args.cells,
        args.genome,
        args.genome_length,
        args.mutation_rate,
        args.steps,
        args.ratio,
    )
