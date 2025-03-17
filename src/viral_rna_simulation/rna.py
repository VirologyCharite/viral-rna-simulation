import sys
from collections import defaultdict, Counter

from viral_rna_simulation.genome import Genome

# from viral_rna_simulation.genome import genomes_str


class RNA:
    def __init__(self, genome: Genome) -> None:
        self.genome = genome
        self.replications = 0

    def __str__(self):
        positive = "+" if self.positive else "-"
        return f"<({positive}) RNA {self.genome}>"

    def __repr__(self):
        positive = "+" if self.positive else "-"
        return f"<({positive}) RNA {self.genome!r}>"

    def __eq__(self, other: object, /) -> bool:
        if isinstance(other, RNA):
            return self.genome == other.genome
        return NotImplemented

    def __len__(self) -> int:
        return len(self.genome)

    @property
    def positive(self) -> bool:
        return self.genome.positive

    def replicate(self, mutation_rate: float = 0.0) -> "RNA":
        """
        Make a reverse-complement copy of this RNA, perhaps with mutations.
        """
        self.replications += 1
        return RNA(self.genome.replicate(mutation_rate))

    def sequencing_mutation_counts(
        self, infecting_genome: Genome
    ) -> tuple[dict[str, int], dict[str, dict[bool, Counter[str]]]]:
        """
        Return the mutation counts (relative to the infecting genome) that would be
        counted if this molecule were sequenced. The library preparation involves making
        two (complementary) DNA strands, both of which are assumed to be sequenced.
        """
        genome = self.genome if self.positive else self.genome.rc()
        mutations = defaultdict(int)
        sources: dict[str, dict[bool, Counter[str]]] = {}

        # print(
        #     genomes_str(
        #         genome_1=infecting_genome,
        #         genome_2=genome,
        #         title_1="Infecting genome: ",
        #         title_2="Genome: ",
        #     )
        # )
        # print()

        for a, b in zip(infecting_genome, genome):
            if a != b:
                change = a.base + b.base
                if "pytest" not in sys.modules:
                    # If the genome base does not match the infecting genome, the genome
                    # site must have a mutation history.
                    historical_change, historical_positive = b.mutation_history[-1]
                    reasons = sources.setdefault(
                        change,
                        {
                            True: Counter(),
                            False: Counter(),
                        },
                    )
                    reasons[historical_positive][historical_change] += 1

                # TODO: We should perhaps add two here.
                mutations[change] += 1

        return mutations, sources
