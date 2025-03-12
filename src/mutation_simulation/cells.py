from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from typing import Iterator
from collections import Counter

from mutation_simulation.cell import Cell
from mutation_simulation.genome import Genome


def replicate_rnas(cell: Cell, steps: int, mutation_rate: float) -> Cell:
    cell.replicate_rnas(steps, mutation_rate)
    return cell


class Cells:
    """
    Maintain a collection of cells, all of which initially contain the same RNA.
    """

    def __init__(self, n_cells: int, infecting_genome: Genome) -> None:
        self.infecting_genome = infecting_genome
        self.cells = [Cell(infecting_genome) for _ in range(n_cells)]

    def __iter__(self) -> Iterator[Cell]:
        return iter(self.cells)

    def __str__(self) -> str:
        result = [f"<Cells with {len(self.cells)} cells>"]
        for i, cell in enumerate(self.cells):
            result.append(f"  {i + 1}: {cell}")

        return "\n".join(result)

    def replicate(
        self, workers: int | None = None, steps: int = 100, mutation_rate: float = 0.0
    ) -> None:
        """
        Replicate (in parallel) each cell for a given number of steps.

        @param workers: The number of concurrent worker processes to allow in the
            process pool.
        @param steps: The number of replication steps each cell should perform.
        @param mutation_rate: The per-base mutation probability.
        """
        result = []

        with ProcessPoolExecutor(max_workers=workers) as executor:
            for cell in executor.map(
                replicate_rnas,
                self.cells,
                repeat(steps),
                repeat(mutation_rate),
            ):
                result.append(cell)

        self.cells = result

    def mutation_counts(self) -> Counter:
        """
        Add up all mutations in all RNA molecules in all cells.
        """
        mutations = Counter()

        for cell in self.cells:
            for rna in cell:
                for site in rna.genome:
                    if mutant := site.mutant():
                        mutations[mutant] += 1

        return mutations

    def rna_count(self) -> int:
        """
        Get the number of all RNA molecules in all cells.
        """
        return sum(len(cell) for cell in self.cells)

    def replication_count(self) -> int:
        """
        Get the number of RNA molecule replications that occurred.
        """
        return sum(rna.replications for cell in self.cells for rna in cell)

    def summary(self) -> str:
        """
        Return a summary of all cells.
        """
        replications = self.replication_count()
        result = []

        result.append(f"RNA molecules: {self.rna_count()}")
        result.append(f"Total RNA molecule replications: {replications}")

        if mutations := self.mutation_counts():
            total = sum(mutations.values())
            rate = total / (replications * len(self.infecting_genome))
            result.append("Actual mutations:")
            result.append(f"  Total: {total}")
            result.append(f"  Rate: {rate:.6f}")
            result.append(
                "  From/to: "
                + ", ".join(
                    f"{mutation}:{count}"
                    for mutation, count in sorted(mutations.items())
                ),
            )
        else:
            result.append("Mutations: None")

        apparent_mutations = Counter()
        for cell in self:
            for rna in cell:
                apparent_mutations += rna.sequence(self.infecting_genome)

        if apparent_mutations:
            result.append("Apparent mutations:")
            total = sum(apparent_mutations.values())
            result.append(f"  Total: {total}")
            result.append(
                "  From/to: "
                + ", ".join(
                    f"{mutation}:{count}"
                    for mutation, count in sorted(apparent_mutations.items())
                ),
            )
        else:
            result.append("Apparent mutations: None")

        return "\n".join(result)
