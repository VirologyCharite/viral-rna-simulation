from random import choice
from typing import Iterator

from viral_rna_simulation.genome import Genome
from viral_rna_simulation.rna import RNA


class Cell:
    def __init__(self, infecting_genome: Genome) -> None:
        assert infecting_genome.positive
        self.infecting_genome = infecting_genome
        self.rnas = [RNA(infecting_genome)]

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
        self,
        steps: int,
        mutate_in: str = "both",
        mutation_rate: float = 0.0,
        ratio: int = 1,
        chooser=choice,
    ) -> None:
        """
        Repeatedly ('steps' times) choose an RNA molecule at random from this cell,
        replicate it (once if it is a (+) RNA or 'ratio' times, if it's (-) RNA)
        according to the given mutation rate, and add the result to the list of RNAs
        in this cell. Note that replicating the RNA genome results in the reverse
        complement sequence being synthesized.

        @param chooser: A function that works like 'random.choice', to be used to choose
            the RNA molecule to replicate at each repetition. This is just used for
            testing, to allow for control over what would otherwise be random.
        """
        for _ in range(steps):
            rna = chooser(self.rnas)
            if rna.genome.positive:
                # Our chosen molecule is positive, so we're about to make a
                # (-) RNA. If we are only mutating positive strands, we must
                # set the mutation rate to zero.
                rate = 0.0 if mutate_in == "positive" else mutation_rate
                self.rnas.append(rna.replicate(rate))
            else:
                rate = 0.0 if mutate_in == "negative" else mutation_rate
                self.rnas.extend(rna.replicate(rate) for _ in range(ratio))
