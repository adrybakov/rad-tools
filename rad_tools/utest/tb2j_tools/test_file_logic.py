import pytest
from rad_tools.tb2j_tools.file_logic import ExchangeModel


def test_empty_filename():
    with pytest.raises(TypeError):
        ExchangeModel(None)


def test_wrong_filename():
    with pytest.raises(FileNotFoundError):
        ExchangeModel("Ah, music. A magic beyond all we do here!")
