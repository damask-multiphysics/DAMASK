import os

import pytest

import damask

def pytest_addoption(parser):
    parser.addoption("--update",
                     action="store_true",
                     default=False)

@pytest.fixture
def update(request):
    """store current results as new reference results."""
    return request.config.getoption("--update")

@pytest.fixture
def reference_dir_base():
    """directory containing reference results."""
    env = damask.Environment()
    return os.path.join(env.rootDir(),'python','tests','reference')
