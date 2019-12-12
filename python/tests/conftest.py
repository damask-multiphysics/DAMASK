import os

import pytest

def pytest_addoption(parser):
    parser.addoption("--update",
                     action="store_true",
                     default=False)

@pytest.fixture
def update(request):
    """Store current results as new reference results."""
    return request.config.getoption("--update")

@pytest.fixture
def reference_dir_base():
    """Directory containing reference results."""
    return os.path.join(os.path.dirname(__file__),'reference')
