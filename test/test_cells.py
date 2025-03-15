from viral_rna_simulation.cells import Cells
from viral_rna_simulation.genome import Genome


class Test_cells:
    """
    Test the Cells class.
    """

    def test_cell_with_one_cell(self) -> None:
        """
        A new Cells with one cell has length one.
        """
        cells = Cells(1, Genome("A"))
        assert len(cells) == 1

    def test_length(self) -> None:
        """
        The Cells __len__ method must work.
        """
        assert len(Cells(5, Genome("ATTC"))) == 5

    def test_infecting_genome_in_five_cells(self) -> None:
        """
        A new Cells should have the expected genome in the RNA for each of its cells.
        """
        cells = Cells(5, Genome("ATTC"))
        genomes = [
            "".join(site.base for site in rna.genome) for cell in cells for rna in cell
        ]
        assert genomes == ["ATTC"] * 5

    def test_cell_with_one_cell_has_one_RNA(self) -> None:
        """
        A new Cells with one cell has one RNA in it.
        """
        cells = Cells(1, Genome("A"))
        assert cells.rna_count() == 1

    def test_no_mutations(self) -> None:
        """
        A new Cells has no mutations.
        """
        cells = Cells(1, Genome("A"))
        assert cells.mutation_counts() == {}

    def test_no_replications(self) -> None:
        """
        A new Cells has not replicated.
        """
        cells = Cells(1, Genome("A"))
        assert cells.replication_count() == 0

    def test_summary(self) -> None:
        """
        Calling the summary method should not result in any error and should return a
        non-empty string.
        """
        cells = Cells(1, Genome("A"))
        summary = cells.summary()
        assert isinstance(summary, str)
        assert summary

    def test_replicate(self) -> None:
        """
        Test replication does not fail.
        """
        cells = Cells(2, Genome("A"))
        assert len(cells) == 2

        cells.replicate(steps=1)
        assert cells.rna_count() == 4
        assert cells.replication_count() == 2

        # At each step, each cell picks one RNA to replicate, so in each call
        # to replicated, the overall number of RNA molecules (i.e., summed
        # over all cells) goes up by the product of the number of workers and
        # the number of steps (2 x 3 = 6, in this call).
        cells.replicate(steps=3)
        assert cells.rna_count() == 10
        assert cells.replication_count() == 8

        cells.replicate(steps=1)
        assert cells.rna_count() == 12
        assert cells.replication_count() == 10

        # No mutations have occurred because the mutation_rate defaults to 0.0
        assert cells.mutation_counts() == {}

        # There should still be two cells.
        assert len(cells) == 2

        # With the expected number of RNA molecules in each.
        assert [len(cell) for cell in cells] == [6, 6]

