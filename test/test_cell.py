import pytest

from viral_rna_simulation.cell import Cell
from viral_rna_simulation.genome import Genome
from viral_rna_simulation.rna import RNA


class Test_basic:
    """
    Test basic properties of the Cell class.
    """

    def test_length(self) -> None:
        """
        A new cell has just one RNA in it.
        """
        assert len(Cell(Genome("A"))) == 1


class Test_replicate:
    """
    Test replication of the Cell class.
    """

    @pytest.mark.parametrize("base", "ACGT")
    def test_count_increases_by_one(self, base) -> None:
        """
        The number of RNA molecules must increase by one after replication.
        """
        cell = Cell(Genome(base))
        assert len(cell) == 1
        cell.replicate_rnas(1)
        assert len(cell) == 2

    @pytest.mark.parametrize("base", "ACGT")
    def test_count_increase_with_ratio(self, base) -> None:
        """
        The number of RNA molecules must increase by one after replication to
        make the negative strand and then by the 'ratio' in making positive
        strands from the negative.
        """
        ratio = 10
        cell = Cell(Genome(base))
        assert len(cell) == 1
        cell.replicate_rnas(1)
        assert len(cell) == 2

        def choose_negative(rnas: list[RNA]):
            infecting, negative = rnas
            assert infecting.positive and not negative.positive
            return negative

        cell.replicate_rnas(1, ratio=ratio, chooser=choose_negative)
        assert len(cell) == ratio + 2
        assert [rna.positive for rna in cell] == [True, False] + [True] * ratio

    @pytest.mark.parametrize("base", "ACGT")
    def test_initial_positive(self, base) -> None:
        """
        The first RNA molecule in a cell must be positive sense.
        """
        cell = Cell(Genome(base))
        (rna,) = cell
        assert rna.positive

    @pytest.mark.parametrize("base", "ACGT")
    def test_replicate_is_negative(self, base) -> None:
        """
        Replication must create a negative RNA.
        """
        cell = Cell(Genome(base))
        cell.replicate_rnas(1)
        _, rna2 = cell
        assert not rna2.positive

    @pytest.mark.parametrize("base", "ACGT")
    def test_replicate_is_reverse_complement(self, base) -> None:
        """
        Replication must create the reverse complement.
        """
        cell = Cell(Genome(base))
        cell.replicate_rnas(1)
        rna1, rna2 = cell
        assert rna1.genome == rna2.genome.rc()
