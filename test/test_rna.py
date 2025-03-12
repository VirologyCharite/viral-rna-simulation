import pytest

from mutation_simulation.genome import Genome
from mutation_simulation.rna import RNA
from mutation_simulation.utils import rc1


class Test_replicate:
    """
    Test the RNA replicate function.
    """

    def test_positive_to_negative(self):
        assert RNA(Genome("A"), True).replicate() == RNA(Genome("T"), False)

    def test_negative_to_positive(self):
        assert RNA(Genome("A"), False).replicate() == RNA(Genome("T"), True)

    def test_roundtrip_positive(self):
        rna = RNA(Genome("GAT"), True)
        assert rna.replicate().replicate() == rna

    def test_roundtrip_negative(self):
        rna = RNA(Genome("GAT"), False)
        assert rna.replicate().replicate() == rna

    def test_three_bases(self):
        assert RNA(Genome("GAT"), True).replicate() == RNA(Genome("ATC"), False)


class Test_equality:
    """
    Tests for equality.
    """

    def test_equals(self):
        rna = RNA(Genome("GAT"), True)
        assert rna == rna

    def test_not_equals_genome(self):
        assert RNA(Genome("CCC"), True) != RNA(Genome("GAT"), True)

    def test_not_equals_parity(self):
        assert RNA(Genome("CCC"), True) != RNA(Genome("CCC"), False)


single_changes = [(a, b) for a in "ACGT" for b in "ACGT"]

# single_changes = [
#     ("A", "C"),
#     ("A", "G"),
#     ("A", "T"),
#     ("C", "G"),
#     ("C", "T"),
#     ("G", "T"),
# ]

# single_changes.extend([(b, a) for (a, b) in single_changes])


class Test_sequence:
    """
    Tests of the 'sequence' method.
    """

    def test_no_changes(self):
        rna = RNA(Genome("GAT"), True)
        mutations = rna.sequencing_counts(Genome("GAT"))
        assert not mutations

    def test_AG_to_CT_positive(self):
        rna = RNA(Genome("CT"), True)
        mutations = rna.sequencing_counts(Genome("AG"))
        expected = {"AC": 1, "GT": 1}
        assert mutations == expected

    @pytest.mark.parametrize("from_,to", single_changes)
    def test_one_change_positive(self, from_, to):
        rna = RNA(Genome(to), True)
        mutations = rna.sequencing_counts(Genome(from_))
        expected = {} if from_ == to else {from_ + to: 1}
        assert mutations == expected

    @pytest.mark.parametrize("from_,to", single_changes)
    def test_two_changes_positive(self, from_, to):
        rna = RNA(Genome(to + to), True)
        mutations = rna.sequencing_counts(Genome(from_ + from_))
        expected = {} if from_ == to else {from_ + to: 2}
        assert mutations == expected

    @pytest.mark.parametrize("from_,to", single_changes)
    def test_one_change_negative(self, from_, to):
        rna = RNA(Genome(to), False)
        mutations = rna.sequencing_counts(Genome(from_))
        expected = {} if from_ == rc1(to) else {from_ + rc1(to): 1}
        assert mutations == expected

    @pytest.mark.parametrize("from_,to", single_changes)
    def test_two_changes_negative(self, from_, to):
        rna = RNA(Genome(to + to), False)
        mutations = rna.sequencing_counts(Genome(from_ + from_))
        expected = {} if from_ == rc1(to) else {from_ + rc1(to): 2}
        assert mutations == expected
