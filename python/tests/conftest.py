from pathlib import Path

import pytest

# Use to monkeypatch damask.version (for comparsion to reference files that contain version information)
def pytest_configure():
    pytest.dummy_version = '99.99.99-9999-pytest'


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
    return Path(__file__).parent/'reference'
