from collections import Counter

from mutation_simulation.genome import Genome


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

    def replicate(self, mutation_rate: float = 0.0) -> "RNA":
        """
        Make a reverse-complement copy of this RNA, perhaps with mutations.
        """
        self.replications += 1
        return RNA(
            self.genome.replicate(mutation_rate),
            not self.positive,
        )

    def orig_sequence(self, infecting_genome: Genome) -> Counter:
        """
        Return the mutations (relative to the infecting genome) that would be counted
        if this molecule were sequenced. The library preparation involves making two
        (complementary) DNA strands, both of which are assumeed to be sequenced.
        """
        counts = []

        print("".join(site.base for site in infecting_genome))
        for genome in self.genome, self.genome.rc():
            print("".join(site.base for site in genome))
            these_counts = Counter()
            for a, b in zip(infecting_genome, genome):
                if a != b:
                    these_counts[a.base + b.base] += 1

            counts.append(these_counts)

        if counts[0] != counts[1]:
            print("Counts are not identical!")
            for these_counts in counts:
                print(
                    f"Total {sum(these_counts.values())}: " +
                    ", ".join(
                        f"{mutation}:{count}"
                        for mutation, count in sorted(these_counts.items())
                    ),
                )

        return counts[0] + counts[1]

    def sequencing_counts(self, infecting_genome: Genome) -> Counter:
        """
        Return the mutation counts (relative to the infecting genome) that would be
        counted if this molecule were sequenced. The library preparation involves making
        two (complementary) DNA strands, both of which are assumed to be sequenced.
        """
        genome = self.genome if self.positive else self.genome.rc()

        mutations = Counter()

        for a, b in zip(infecting_genome, genome):
            if a != b:
                # TODO: We should add two here.
                mutations[a.base + b.base] += 1

        return mutations
