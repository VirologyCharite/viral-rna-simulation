from viral_rna_simulation.genome import Genome


class Test_basic:
    """
    Test basic properties of the Genome class.
    """

    def test_length(self) -> None:
        """
        A genome must have the expected length.
        """
        assert len(Genome("AG")) == 2

    def test_str(self) -> None:
        """
        The str method must work.
        """
        assert str(Genome("AG")) == "AG"

    def test_equality(self) -> None:
        """
        The __eq__ method must work when genomes are equal.
        """
        assert Genome("AG") == Genome("AG")

    def test_inequality(self) -> None:
        """
        The __eq__ method must work when genomes are not equal.
        """
        assert Genome("AG") != Genome("GA")

    def test_equality_different_polarity(self) -> None:
        """
        The __eq__ method must work when genomes nucleotides are equal but the genomes
        are of different polarity.
        """
        assert Genome("AG") != Genome("AG", positive=False)

    def test_getitem(self) -> None:
        """
        The __getitem__ method must work.
        """
        print(repr(Genome("AG")))
        assert Genome("AG")[0].base == "A"
        assert Genome("AG")[1].base == "G"
