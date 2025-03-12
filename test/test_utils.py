import pytest

from mutation_simulation.utils import rc


class Test_rc:
    """
    Test the rc function.
    """
    @pytest.mark.parametrize(
        "from_,to",
        [
            ("A", "T"),
            ("T", "A"),
            ("C", "G"),
            ("G", "C"),
        ],
    )
    def test_one_char(self, from_, to):
        assert rc(from_) == to

    @pytest.mark.parametrize(
        "from_,to",
        [
            ("ATCG", "CGAT"),
            ("AAC", "GTT"),
        ],
    )
    def test_several_chars(self, from_, to):
        assert rc(from_) == to
