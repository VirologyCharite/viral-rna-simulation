import pytest

from viral_rna_simulation.utils import rc, rc1, mutate_base


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
    def test_one_char(self, from_, to) -> None:
        assert rc(from_) == to

    @pytest.mark.parametrize(
        "from_,to",
        [
            ("ATCG", "CGAT"),
            ("AAC", "GTT"),
        ],
    )
    def test_several_chars(self, from_, to) -> None:
        assert rc(from_) == to


class Test_rc1:
    """
    Test the rc1 function.
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
    def test_one_base(self, from_, to) -> None:
        assert rc1(from_) == to


class Test_mutate_base:
    @pytest.mark.parametrize(
        "from_,expected",
        [
            ("A", "CGT"),
            ("C", "AGT"),
            ("G", "ACT"),
            ("T", "ACG"),
        ],
    )
    def test_one_base(self, from_, expected) -> None:
        for _ in range(100):
            assert mutate_base(from_) in expected
