from functools import cache
from random import choice


COMPLEMENT = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
}

MUTANTS = {
    "A": "CGT",
    "C": "AGT",
    "G": "ACT",
    "T": "ACG",
}


@cache
def rc(s: str) -> str:
    return "".join(COMPLEMENT[base] for base in reversed(s))

def rc1(s: str) -> str:
    return COMPLEMENT[s]


def mutate_base(base: str) -> str:
    return choice(MUTANTS[base])
