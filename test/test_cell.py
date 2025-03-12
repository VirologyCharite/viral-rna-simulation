import pytest

from mutation_simulation.genome import Genome
from mutation_simulation.cell import Cell


class Test_basic:
    """
    Test basic properties of the Cell class.
    """

    def test_length(self):
        """
        A new cell has just one RNA in it.
        """
        assert len(Cell(Genome("A"))) == 1


class Test_replicate:
    """
    Test replication of the Cell class.
    """

    def test_length(self):
        """
        The number of molecules must increase by one after replication.
        """
        cell = Cell(Genome("A"))
        assert len(cell) == 1
        cell.replicate_rnas(1)
        assert len(cell) == 2

    def test_initial_positive(self):
        """
        The first RNA molecule in a cell must be positive sense.
        """
        cell = Cell(Genome("A"))
        (rna,) = cell
        assert rna.positive

    def test_replicate_is_negative(self):
        """
        Replication must create a negative RNA.
        """
        cell = Cell(Genome("A"))
        cell.replicate_rnas(1)
        _, rna2 = cell
        assert not rna2.positive

    def test_replicate_is_reverse_complement(self):
        """
        Replication must create the reverse complement.
        """
        cell = Cell(Genome("A"))
        cell.replicate_rnas(1)
        rna1, rna2 = cell
        assert rna1.genome == rna2.genome.rc()

    def test_second_replicate_is_same_as_original(self):
        """
        Two rounds of replication must recreate the original RNA.
        """
        cell = Cell(Genome("A"))
        cell.replicate_rnas(2)
        rna1, _, rna3 = cell
        assert rna1 == rna3
