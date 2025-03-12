from random import choice
from typing import Iterator

from mutation_simulation.site import Site


class Genome:
    def __init__(self, sites: list[Site] | str | None = None, length: int = 0) -> None:
        if sites:
            self.sites = (
                sites if isinstance(sites, list) else [Site(base) for base in sites]
            )
        elif length:
            self.sites = [Site(choice("ACGT")) for _ in range(length)]
        else:
            raise ValueError(
                "You must provide either the genome sites or a non-zero genome length."
            )

    def __iter__(self) -> Iterator[Site]:
        return iter(self.sites)

    def __len__(self) -> int:
        return len(self.sites)

    def __eq__(self, other: object, /) -> bool:
        if isinstance(other, Genome):
            return all(a == b for (a, b) in zip(self.sites, other.sites))
        return NotImplemented

    def __str__(self) -> str:
        result = [f"<Genome length {len(self.sites)}>"]

        for site in self.sites:
            result.append(f"  {site}")

        return "\n".join(result)

    def replicate(self, mutation_rate: float = 0.0) -> "Genome":
        """
        Copy the new genome (reverse complemented), possibly with mutations.
        """
        return Genome([site.replicate(mutation_rate) for site in reversed(self.sites)])

    def rc(self) -> "Genome":
        """
        Get a reverse-complemented genome.
        """
        return Genome([site.rc() for site in reversed(self.sites)])
