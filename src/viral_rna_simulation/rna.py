from collections import defaultdict

from viral_rna_simulation.genome import Genome

# from viral_rna_simulation.genome import genomes_str


class RNA:
    def __init__(
        self,
        genome: Genome,
        positive: bool,
    ) -> None:
        self.genome = genome
        self.positive = positive
        self.replications = 0

    def __str__(self):
        polarity = "+" if self.positive else "-"
        return f"<{polarity}RNA {self.genome}>"

    def __eq__(self, other: object, /) -> bool:
        if isinstance(other, RNA):
            return self.genome == other.genome and self.positive == other.positive
        return NotImplemented

    def __len__(self) -> int:
        return len(self.genome)

    def replicate(self, mutation_rate: float = 0.0) -> "RNA":
        """
        Make a reverse-complement copy of this RNA, perhaps with mutations.
        """
        self.replications += 1
        return RNA(
            self.genome.replicate(mutation_rate),
            not self.positive,
        )

    def sequencing_mutation_counts(self, infecting_genome: Genome) -> dict[str, int]:
        """
        Return the mutation counts (relative to the infecting genome) that would be
        counted if this molecule were sequenced. The library preparation involves making
        two (complementary) DNA strands, both of which are assumed to be sequenced.
        """
        genome = self.genome if self.positive else self.genome.rc()
        mutations = defaultdict(int)

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
                # TODO: We should perhaps add two here.
                mutations[a.base + b.base] += 1

        return mutations
