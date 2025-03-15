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
        assert RNA(Genome("A"), True).replicate() == RNA(Genome("T"), False)

    def test_negative_to_positive(self) -> None:
        assert RNA(Genome("A"), False).replicate() == RNA(Genome("T"), True)

    def test_roundtrip_positive(self) -> None:
        rna = RNA(Genome("GAT"), True)
        assert rna.replicate().replicate() == rna

    def test_roundtrip_negative(self) -> None:
        rna = RNA(Genome("GAT"), False)
        assert rna.replicate().replicate() == rna

    def test_three_bases(self) -> None:
        assert RNA(Genome("GAT"), True).replicate() == RNA(Genome("ATC"), False)


class Test_equality:
    """
    Tests for equality.
    """

    def test_equals(self) -> None:
        rna = RNA(Genome("GAT"), True)
        assert rna == rna

    def test_not_equals_genome(self) -> None:
        assert RNA(Genome("CCC"), True) != RNA(Genome("GAT"), True)

    def test_not_equals_parity(self) -> None:
        assert RNA(Genome("CCC"), True) != RNA(Genome("CCC"), False)


single_changes = [(a, b) for a in "ACGT" for b in "ACGT"]


class Test_sequencing_mutation_counts:
    """
    Tests of the 'sequencing_mutation_countssequence' method.
    """

    def test_no_changes(self) -> None:
        rna = RNA(Genome("GAT"), True)
        mutations = rna.sequencing_mutation_counts(Genome("GAT"))
        assert not mutations

    def test_AG_to_CT_positive(self) -> None:
        rna = RNA(Genome("CT"), True)
        mutations = rna.sequencing_mutation_counts(Genome("AG"))
        expected = {"AC": 1, "GT": 1}
        assert mutations == expected

    @pytest.mark.parametrize("from_,to", single_changes)
    def test_one_change_positive(self, from_, to) -> None:
        rna = RNA(Genome(to), True)
        mutations = rna.sequencing_mutation_counts(Genome(from_))
        expected = {} if from_ == to else {from_ + to: 1}
        assert mutations == expected

    @pytest.mark.parametrize("from_,to", single_changes)
    def test_two_changes_positive(self, from_, to) -> None:
        rna = RNA(Genome(to + to), True)
        mutations = rna.sequencing_mutation_counts(Genome(from_ + from_))
        expected = {} if from_ == to else {from_ + to: 2}
        assert mutations == expected

    @pytest.mark.parametrize("from_,to", single_changes)
    def test_one_change_negative(self, from_, to) -> None:
        rna = RNA(Genome(to), False)
        mutations = rna.sequencing_mutation_counts(Genome(from_))
        expected = {} if from_ == rc1(to) else {from_ + rc1(to): 1}
        assert mutations == expected

    @pytest.mark.parametrize("from_,to", single_changes)
    def test_two_changes_negative(self, from_, to) -> None:
        rna = RNA(Genome(to + to), False)
        mutations = rna.sequencing_mutation_counts(Genome(from_ + from_))
        expected = {} if from_ == rc1(to) else {from_ + rc1(to): 2}
        assert mutations == expected

    def test_longer_positive(self) -> None:
        rna = RNA(Genome("AA"), True)
        mutations = rna.sequencing_mutation_counts(Genome("CC"))
        assert mutations == {"CA": 2}

    def test_longer_negative(self) -> None:
        rna = RNA(Genome("AA"), False)
        mutations = rna.sequencing_mutation_counts(Genome("CC"))
        assert mutations == {"CT": 2}

    def test_peter_email_example_1(self) -> None:
        rna = RNA(Genome("A"), False)
        mutations = rna.sequencing_mutation_counts(Genome("G"))
        assert mutations == {"GT": 1}

    def test_peter_email_example_2(self) -> None:
        rna = RNA(Genome("A"), True)
        mutations = rna.sequencing_mutation_counts(Genome("C"))
        assert mutations == {"CA": 1}

    def test_mutation_in_making_the_negative_which_is_then_copied_many_times(
        self,
    ) -> None:
        infecting_genome = Genome("G")
        sequencing_mutations = Counter()

        # Make a negative RNA with an 'A' which is a mutation, since the infecting
        # genome has a 'G' an error-free negative rc copy would have a 'C'.
        negative = RNA(Genome("A"), False)

        sequencing_mutations.update(
            negative.sequencing_mutation_counts(infecting_genome)
        )

        assert sequencing_mutations == {"GT": 1}

        # Copy the negative 10 times (with no error). This will create 10
        # positive RNAs with a 'T' genome.
        positives = [negative.replicate() for _ in range(10)]
        assert all(
            len(rna) == 1 and rna.positive and rna.genome[0].base == "T"
            for rna in positives
        )

        for rna in positives:
            sequencing_mutations.update(
                rna.sequencing_mutation_counts(infecting_genome)
            )

        assert sequencing_mutations == {"GT": 11}

    def test_mutation_in_making_the_positive(self) -> None:
        infecting_genome = Genome("C")
        sequencing_mutations = Counter()

        # The negative, with no mutation.
        negative = RNA(Genome("G"), False)

        sequencing_mutations.update(
            negative.sequencing_mutation_counts(infecting_genome)
        )

        # Make a positive RNA with an 'A' which is a mutation, since the infecting
        # genome has a 'C', an error-free negative rc copy would have a 'G', and so
        # an rc copy of that to make another positive would bring us back to 'C'.
        positive = RNA(Genome("A"), True)

        sequencing_mutations.update(
            positive.sequencing_mutation_counts(infecting_genome)
        )

        assert sequencing_mutations == {"CA": 1}
