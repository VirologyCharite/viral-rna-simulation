from random import uniform

from viral_rna_simulation.utils import mutate_base, rc1


class Site:
    """
    Manage a genome site.

    @param base: The nucleotide.
    @param offset: The 0-based offset of this site within its genome.
    @param mutant: True if the site was created by a mutation (in which case
        the detail of the mutation is in mutation_history[-1]).
    @param mutation_history: A list of 2-tuples that have a nucleotide from/to change
        (e.g, "AT"), followed by a bool that is True if the change took place in a
        positive RNA, else False.
    """
    def __init__(
        self,
        base: str,
        offset: int,
        mutant: bool = False,
        mutation_history: list[tuple[str, bool]] | None = None,
    ) -> None:
        self.base = base
        self.offset = offset
        self.mutant = mutant
        self.mutation_history = mutation_history or []

    def __str__(self) -> str:
        mutations = (
            (
                " Mutations="
                + " ".join(
                    f"{'+' if positive else '-'}{change}"
                    for (change, positive) in self.mutation_history
                )
            )
            if self.mutation_history
            else ""
        )
        return f"<Site offset={self.offset}, mutant={self.mutant}, base={self.base!r}{mutations}>"

    def __eq__(self, other) -> bool:
        if isinstance(other, Site):
            # We should never be comparing sites at different offsets.
            assert self.offset == other.offset
            return self.base == other.base
        return NotImplemented

    def replicate(self, positive: bool, mutation_rate: float = 0.0) -> "Site":
        """
        Make a replicate (in reverse complement) of this site.

        @param positive: The (+/-) state of the new site.
        @param mutation_rate: The mutation rate used to decide whether the new
            site should be a mutant.
        """
        rc_base = rc1(self.base)
        if mutation_rate > 0.0 and uniform(0.0, 1.0) <= mutation_rate:
            mutant = True
            new_base = mutate_base(rc_base)
            change = rc_base + new_base
            # Or: change = self.base + new_base (depends on what we're saying changed).
            mutation_history = self.mutation_history + [(change, positive)]
        else:
            mutant = False
            new_base = rc_base
            mutation_history = self.mutation_history[:]

        return Site(
            new_base, self.offset, mutant=mutant, mutation_history=mutation_history
        )

    def rc(self) -> "Site":
        """
        Return a reverse-complemented site.
        """
        return Site(
            rc1(self.base),
            offset=self.offset,
            mutant=False,
            mutation_history=self.mutation_history[:],
        )
