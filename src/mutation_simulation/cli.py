import argparse

from mutation_simulation.simulate import run


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells", type=int, default=1)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--genome-length", type=int)
    group.add_argument("--genome")
    parser.add_argument("--steps", type=int, default=1000)
    parser.add_argument("--mutation-rate", type=float, default=0.001)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run(args.cells, args.genome, args.genome_length, args.mutation_rate, args.steps)
