from mutation_simulation.cells import Cells
from mutation_simulation.genome import Genome


def run(
    n_cells: int,
    genome: str | None,
    genome_length: int,
    mutation_rate: float,
    steps: int,
) -> None:
    """
    Simulate a number of cells.
    """
    infecting_genome = Genome(genome, genome_length)
    cells = Cells(n_cells, infecting_genome)

    cells.replicate(steps=steps, mutation_rate=mutation_rate)

    print(cells.summary())
