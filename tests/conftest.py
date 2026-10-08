import pytest
import torch

from plm_interp.models.loader import load_model

UBIQUITIN = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"


@pytest.fixture(scope="session")
def bundle():
    torch.set_num_threads(2)
    return load_model("8M")


@pytest.fixture(scope="session")
def seq():
    return UBIQUITIN
