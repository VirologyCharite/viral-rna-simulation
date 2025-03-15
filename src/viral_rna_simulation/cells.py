from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from typing import Iterator
from collections import Counter

from viral_rna_simulation.cell import Cell
from viral_rna_simulation.genome import Genome
from viral_rna_simulation.utils import mutations_str


def replicate_rnas(cell: Cell, steps: int, mutation_rate: float, ratio: int) -> Cell:
    cell.replicate_rnas(steps, mutation_rate=mutation_rate, ratio=ratio)
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

    def __len__(self) -> int:
        return len(self.cells)

    def __str__(self) -> str:
        result = [f"<Cells with {len(self.cells)} cells>"]
        for i, cell in enumerate(self.cells):
            result.append(f"  {i + 1}: {cell}")

        return "\n".join(result)

    def replicate(
            self, workers: int | None = None, steps: int = 1, mutation_rate: float = 0.0, ratio: int = 1
    ) -> None:
        """
        Replicate (in parallel) each cell for a given number of steps.

        At each step, each cell picks one RNA to replicate, so in each call
        to replicated, the number of RNA molecules overall (i.e., summed over
        all cells) goes up by the product of the number of workers and the
        number of steps (2 x 3 = 6, in this call).

        @param workers: The number of concurrent worker processes to allow in the
            process pool.
        @param steps: The number of replication steps each cell should perform.
        @param mutation_rate: The per-base mutation probability.
        @param ratio: The number of +RNA molecules to make from a -RNA.
        """
        new_cells = []

        with ProcessPoolExecutor(max_workers=workers) as executor:
            for cell in executor.map(
                replicate_rnas,
                self.cells,
                repeat(steps),
                repeat(mutation_rate),
                repeat(ratio),
            ):
                new_cells.append(cell)

        self.cells = new_cells

    def mutation_counts(self) -> Counter:
        """
        Add up all mutations in all RNA molecules in all cells.
        """
        mutations = Counter()

        for cell in self.cells:
            for rna in cell:
                for site in rna.genome:
                    mutations.update(site.mutations)

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
            total_mutations = sum(mutations.values())
            result.append("Actual mutations:")
            result.append(f"  Total: {total_mutations}")

            positive_mutation_count = negative_mutation_count = 0
            for cell in self:
                for rna in cell:
                    count = sum(rna.genome.mutations().values())
                    if rna.positive:
                        positive_mutation_count += count
                    else:
                        negative_mutation_count += count

            assert total_mutations == positive_mutation_count + negative_mutation_count
            result.append(
                f"    Of which, {positive_mutation_count} in (+) RNA molecules"
            )
            result.append(
                f"    Of which, {negative_mutation_count} in (-) RNA molecules"
            )

            rate = total_mutations / (replications * len(self.infecting_genome))
            result.append(f"  Rate: {rate:.6f}")
            result.append(
                "  From/to: "
                + ", ".join(
                    f"{mutation}:{count}"
                    for mutation, count in sorted(mutations.items())
                ),
            )
        else:
            total_mutations = 0
            result.append("Mutations: None")

        apparent_mutations = Counter()
        for cell in self:
            for rna in cell:
                apparent_mutations += rna.sequencing_mutation_counts(
                    self.infecting_genome
                )

        if apparent_mutations:
            total = sum(apparent_mutations.values())
            result.extend([
                "Apparent mutations:",
                f"  Total: {total}",
                f"  From/to: {mutations_str(apparent_mutations)}",
            ])
        else:
            result.append("Apparent mutations: None")

        return "\n".join(result)
