from viral_rna_simulation.cells import Cells
from viral_rna_simulation.genome import Genome


def run(
    n_cells: int,
    genome: str | None,
    genome_length: int,
    mutate_in: str,
    mutation_rate: float,
    steps: int,
    ratio: int,
) -> Cells:
    """
    Simulate a number of cells.
    """
    infecting_genome = Genome(genome, genome_length)
    cells = Cells(n_cells, infecting_genome)

    cells.replicate(
        steps=steps, mutate_in=mutate_in, mutation_rate=mutation_rate, ratio=ratio
    )

    return cells
