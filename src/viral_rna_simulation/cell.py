from random import choice
from typing import Iterator

from viral_rna_simulation.genome import Genome
from viral_rna_simulation.rna import RNA


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

    def replicate_rnas(
        self, steps: int, mutation_rate: float = 0.0, ratio: int = 1, chooser=choice
    ) -> None:
        """
        Repeatedly ('steps' times) randomly choose an RNA molecule from this cell,
        replicate it ('ratio' times, if the RNA is a negative strand) according to some
        mutation rate, and add the result to the list of RNAs in this cell. Note that
        replicating the RNA genome results in the reverse complement sequence being
        synthesized.

        @param chooser: A function that works like 'random.choice', to be used to choose
            the RNA molecule to replicate at each repetition. This is just used for
            testing, to allow for control over what would otherwise be random.
        """
        for _ in range(steps):
            rna = chooser(self.rnas)
            if rna.positive:
                self.rnas.append(rna.replicate(mutation_rate))
            else:
                self.rnas.extend(rna.replicate(mutation_rate) for _ in range(ratio))
