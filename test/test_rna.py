import pytest
from collections import Counter

from viral_rna_simulation.genome import Genome
from viral_rna_simulation.rna import RNA
from viral_rna_simulation.utils import rc1


class Test_replicate:
    """
    Test the RNA replicate function.
    """

    def test_positive_to_negative(self) -> None:
        assert RNA(Genome("A", positive=True)).replicate() == RNA(Genome("T", positive=False))

    def test_negative_to_positive(self) -> None:
        assert RNA(Genome("A", positive=False)).replicate() == RNA(Genome("T", positive=True))

    def test_roundtrip_positive(self) -> None:
        rna = RNA(Genome("GAT", positive=True))
        assert rna.replicate().replicate() == rna

    def test_roundtrip_negative(self) -> None:
        rna = RNA(Genome("GAT", positive=False))
        assert rna.replicate().replicate() == rna

    def test_three_bases(self) -> None:
        assert RNA(Genome("GAT", positive=True)).replicate() == RNA(Genome("ATC", positive=False))


class Test_equality:
    """
    Tests for equality.
    """

    def test_equals(self) -> None:
        rna = RNA(Genome("GAT", positive=True))
        assert rna == rna

    def test_not_equals_genome(self) -> None:
        assert RNA(Genome("CCC", positive=True)) != RNA(Genome("GAT", positive=True))

    def test_not_equals_parity(self) -> None:
        assert RNA(Genome("CCC", positive=True)) != RNA(Genome("CCC", positive=False))


single_changes = [(a, b) for a in "ACGT" for b in "ACGT"]


class Test_sequencing_mutation_counts:
    """
    Tests of the 'sequencing_mutation_countssequence' method.
    """

    def test_no_changes(self) -> None:
        rna = RNA(Genome("GAT", positive=True))
        mutations, _ = rna.sequencing_mutation_counts(Genome("GAT"))
        assert not mutations

    def test_AG_to_CT_positive(self) -> None:
        rna = RNA(Genome("CT", positive=True))
        mutations, _ = rna.sequencing_mutation_counts(Genome("AG"))
        expected = {"AC": 1, "GT": 1}
        assert mutations == expected

    @pytest.mark.parametrize("from_,to", single_changes)
    def test_one_change_positive(self, from_, to) -> None:
        rna = RNA(Genome(to, True))
        mutations, _ = rna.sequencing_mutation_counts(Genome(from_))
        expected = {} if from_ == to else {from_ + to: 1}
        assert mutations == expected

    @pytest.mark.parametrize("from_,to", single_changes)
    def test_two_changes_positive(self, from_, to) -> None:
        rna = RNA(Genome(to + to, True))
        mutations, _ = rna.sequencing_mutation_counts(Genome(from_ + from_))
        expected = {} if from_ == to else {from_ + to: 2}
        assert mutations == expected

    @pytest.mark.parametrize("from_,to", single_changes)
    def test_one_change_negative(self, from_, to) -> None:
        rna = RNA(Genome(to, positive=False))
        mutations, _ = rna.sequencing_mutation_counts(Genome(from_))
        expected = {} if from_ == rc1(to) else {from_ + rc1(to): 1}
        assert mutations == expected

    @pytest.mark.parametrize("from_,to", single_changes)
    def test_two_changes_negative(self, from_, to) -> None:
        rna = RNA(Genome(to + to, positive=False))
        mutations, _ = rna.sequencing_mutation_counts(Genome(from_ + from_))
        expected = {} if from_ == rc1(to) else {from_ + rc1(to): 2}
        assert mutations == expected

    def test_longer_positive(self) -> None:
        rna = RNA(Genome("AA", positive=True))
        mutations, _ = rna.sequencing_mutation_counts(Genome("CC"))
        assert mutations == {"CA": 2}

    def test_longer_negative(self) -> None:
        rna = RNA(Genome("AA", positive=False))
        mutations, _ = rna.sequencing_mutation_counts(Genome("CC"))
        assert mutations == {"CT": 2}

    def test_peter_email_example_1(self) -> None:
        rna = RNA(Genome("A", positive=False))
        mutations, _ = rna.sequencing_mutation_counts(Genome("G"))
        assert mutations == {"GT": 1}

    def test_peter_email_example_2(self) -> None:
        rna = RNA(Genome("A", positive=True))
        mutations, _ = rna.sequencing_mutation_counts(Genome("C"))
        assert mutations == {"CA": 1}

    def test_mutation_in_making_the_negative_which_is_then_copied_many_times(
        self,
    ) -> None:
        infecting_genome = Genome("G")
        sequencing_mutations = Counter()

        # Make a negative RNA with an 'A' which is a mutation, since the infecting
        # genome has a 'G' an error-free negative rc copy would have a 'C'.
        negative = RNA(Genome("A", positive=False))

        mutations, _ = negative.sequencing_mutation_counts(infecting_genome)
        sequencing_mutations.update(mutations)

        assert sequencing_mutations == {"GT": 1}

        # Copy the negative 10 times (with no error). This will create 10
        # positive RNAs with a 'T' genome.
        positives = [negative.replicate() for _ in range(10)]
        assert all(
            len(rna) == 1 and rna.positive and rna.genome[0].base == "T"
            for rna in positives
        )

        for rna in positives:
            mutations, _ = rna.sequencing_mutation_counts(infecting_genome)
            sequencing_mutations.update(mutations)

        assert sequencing_mutations == {"GT": 11}

    def test_mutation_in_making_the_positive(self) -> None:
        infecting_genome = Genome("C")
        sequencing_mutations = Counter()

        # The negative, with no mutation.
        negative = RNA(Genome("G", positive=False))

        mutations, _ = negative.sequencing_mutation_counts(infecting_genome)
        sequencing_mutations.update(mutations)

        # Make a positive RNA with an 'A' which is a mutation, since the infecting
        # genome has a 'C', an error-free negative rc copy would have a 'G', and so
        # an rc copy of that to make another positive would bring us back to 'C'.
        positive = RNA(Genome("A", positive=True))

        mutations, _ = positive.sequencing_mutation_counts(infecting_genome)
        sequencing_mutations.update(mutations)

        assert sequencing_mutations == {"CA": 1}
