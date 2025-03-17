from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from typing import Iterator
from collections import Counter

from viral_rna_simulation.cell import Cell
from viral_rna_simulation.genome import Genome
from viral_rna_simulation.utils import mutations_str


def replicate_rnas(
    cell: Cell, steps: int, mutate_in: str, mutation_rate: float, ratio: int
) -> Cell:
    cell.replicate_rnas(
        steps, mutate_in=mutate_in, mutation_rate=mutation_rate, ratio=ratio
    )
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
        self,
        workers: int | None = None,
        steps: int = 1,
        mutate_in: str = "both",
        mutation_rate: float = 0.0,
        ratio: int = 1,
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
        @param mutate_in: The type of RNA molecules to allow mutations in. If
            'negative' or 'positive', mutations should only be allowed in those
            molecules.
        @param mutation_rate: The per-base mutation probability.
        @param ratio: The number of +RNA molecules to make from a -RNA.
        """
        new_cells = []

        with ProcessPoolExecutor(max_workers=workers) as executor:
            for cell in executor.map(
                replicate_rnas,
                self.cells,
                repeat(steps),
                repeat(mutate_in),
                repeat(mutation_rate),
                repeat(ratio),
            ):
                new_cells.append(cell)

        self.cells = new_cells

    def mutation_counts(self) -> tuple[Counter, Counter]:
        """
        Add up all mutations in all (+/-) RNA molecules in all cells.
        """
        positive_counts = Counter()
        negative_counts = Counter()

        for cell in self.cells:
            for rna in cell:
                for site in rna.genome:
                    if site.mutant:
                        change, positive = site.mutation_history[-1]
                        mutations = positive_counts if positive else negative_counts
                        mutations[change] += 1

        return positive_counts, negative_counts

    def rna_count(self) -> tuple[int, int]:
        """
        Get the number of all (+/-) RNA molecules in all cells.
        """
        positive = negative = 0

        for cell in self.cells:
            for rna in cell:
                if rna.genome.positive:
                    positive += 1
                else:
                    negative += 1

        return positive, negative

    def replication_count(self) -> tuple[int, int]:
        """
        Get the number of (+/-) RNA molecule replications that occurred.
        """
        positive = negative = 0

        for cell in self.cells:
            for rna in cell:
                if rna.genome.positive:
                    positive += rna.replications
                else:
                    negative += rna.replications

        return positive, negative

    def apparent_mutation_counts(self) -> tuple[Counter, Counter]:
        """
        Get the apparent changes. I.e., what it looks like happened, based on sample
        preparation, sequencing, alignment to the (+) RNA reference (infecting) genome.
        """
        from_positive = Counter()
        from_negative = Counter()

        for cell in self:
            for rna in cell:
                mutations = from_positive if rna.genome.positive else from_negative
                counts, reasons = rna.sequencing_mutation_counts(self.infecting_genome)
                mutations += counts

        return from_positive, from_negative

    def summary(self) -> str:
        """
        Return a summary of all cells for printing.
        """
        result = []

        # RNA counts.
        positive_rna_count, negative_rna_count = self.rna_count()
        overall_rna_count = positive_rna_count + negative_rna_count
        result.append(f"RNA molecules: {overall_rna_count}")
        result.append(f"  (+) {positive_rna_count}")
        result.append(f"  (-) {negative_rna_count}")

        # Replication counts.
        positive_replications, negative_replications = self.replication_count()
        overall_replications = positive_replications + negative_replications
        result.append(f"Total RNA molecule replications: {overall_replications}")
        result.append(f"  (+): {positive_replications}")
        result.append(f"  (-): {negative_replications}")

        # Changes.
        positive_changes, negative_changes = self.mutation_counts()
        overall_changes = positive_changes + negative_changes

        if overall_changes:
            length = len(self.infecting_genome)

            total_change_count = sum(overall_changes.values())
            positive_change_count = sum(positive_changes.values())
            negative_change_count = sum(negative_changes.values())

            overall_rate = total_change_count / (overall_replications * length)
            positive_rate = positive_change_count / (positive_replications * length)
            negative_rate = negative_change_count / (negative_replications * length)

            result.extend([
                f"Mutations: {total_change_count}",
                f"  In (+) RNA: {positive_change_count}",
                f"  In (-) RNA: {negative_change_count}",
                f"  Overall rate: {overall_rate:.6f}",
                f"    (+) RNA: {positive_rate:.6f}",
                f"    (-) RNA: {negative_rate:.6f}",
                f"  From/to: {mutations_str(overall_changes)}",
            ])

            if positive_changes:
                result.append(f"    (+) RNA: {mutations_str(positive_changes)}")
            if negative_changes:
                result.append(f"    (-) RNA: {mutations_str(negative_changes)}")
        else:
            result.append("Mutations: None")

        from_positive, from_negative = self.apparent_mutation_counts()
        apparent_changes = from_positive + from_negative
        if apparent_changes:
            total = sum(apparent_changes.values())
            result.extend([
                "Apparent mutations:",
                f"  Total: {total}",
                f"  From/to: {mutations_str(apparent_changes)}",
            ])
            if from_positive:
                result.append(f"    (+) From/to: {mutations_str(from_positive)}")
            if from_negative:
                result.append(f"    (-) From/to: {mutations_str(from_negative)}")
        else:
            result.append("Apparent mutations: None")

        return "\n".join(result)
