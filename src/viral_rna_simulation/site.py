from random import uniform
from collections import defaultdict

from viral_rna_simulation.utils import mutate_base, rc1


class Site:
    def __init__(
        self,
        base: str,
        mutations: dict[str, int] | None = None,
    ) -> None:
        self.base = base
        self.mutations = mutations or defaultdict(int)

    def __str__(self) -> str:
        if self.mutations:
            mutations = " " + ", ".join(
                f"{mutation}:{count}"
                for mutation, count in sorted(self.mutations.items())
            )
        else:
            mutations = ""
        return f"<Site {self.base!r}{mutations}>"

    def __eq__(self, other) -> bool:
        if isinstance(other, Site):
            return self.base == other.base
        return NotImplemented

    def replicate(self, mutation_rate: float = 0.0) -> "Site":
        rc_base = rc1(self.base)
        mutations = defaultdict(int)
        if mutation_rate > 0.0 and uniform(0.0, 1.0) <= mutation_rate:
            new_base = mutate_base(rc_base)
            mutations[rc_base + new_base] = 1
            # mutations[self.base + new_base] = 1
        else:
            new_base = rc_base

        return Site(new_base, mutations=mutations)

    def rc(self) -> "Site":
        """
        Return a reverse-complemented site.
        """
        return Site(rc1(self.base))
