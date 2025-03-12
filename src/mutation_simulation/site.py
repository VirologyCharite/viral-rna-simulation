from random import uniform

from mutation_simulation.utils import mutate_base, rc1


class Site:
    def __init__(self, base: str, original_base: str | None = None) -> None:
        self.base = base
        self.original_base = original_base or base

    def __str__(self) -> str:
        was = "" if self.original_base == self.base else f" (was {self.original_base})"
        return f"<Site {self.base!r}{was}>"

    def __eq__(self, other) -> bool:
        if isinstance(other, Site):
            return self.base == other.base
        return NotImplemented

    def mutant(self) -> str | None:
        """
        Is this site a mutant? If so, return the mutation.
        """
        return (
            None if self.original_base == self.base else self.original_base + self.base
        )

    def replicate(self, mutation_rate: float = 0.0) -> "Site":
        rc_base = rc1(self.base)
        if mutation_rate > 0.0 and uniform(0.0, 1.0) <= mutation_rate:
            new_base = mutate_base(rc_base)
        else:
            new_base = rc_base

        return Site(new_base, original_base=rc_base)

    def rc(self) -> "Site":
        """
        Return a reverse-complemented site.
        """
        return Site(rc1(self.base))
