from random import choice
from typing import Iterator

from mutation_simulation.genome import Genome
from mutation_simulation.rna import RNA


class Cell:
    def __init__(self, infecting_genome: Genome) -> None:
        self.infecting_genome = infecting_genome
        self.rnas = [RNA(infecting_genome, positive=True)]

    def __iter__(self) -> Iterator[RNA]:
        return iter(self.rnas)

    def __len__(self) -> int:
        return len(self.rnas)

    def __str__(self) -> str:
        result = [f"<Cell with {len(self.rnas)} RNA molecules>"]
        for i, rna in enumerate(self.rnas):
            result.append(f"    {i + 1}: {rna}")

        return "\n".join(result)

    def replicate_rnas(self, steps: int, mutation_rate: float = 0.0) -> None:
        for _ in range(steps):
            rna = choice(self.rnas)
            self.rnas.append(rna.replicate(mutation_rate))
