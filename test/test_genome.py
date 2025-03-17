import pytest

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
