# Re-export the primary API and formatting helpers from the package
# using relative imports so ``import faststat`` works correctly.
from .faststat import *
from .format import *
from .profiler import Profiler, cli as profiler_cli

try:
    from importlib.metadata import version as _version, PackageNotFoundError
except ImportError:  # pragma: no cover - Python <3.8
    from importlib_metadata import version as _version, PackageNotFoundError

try:
    __version__ = _version(__name__)
except PackageNotFoundError:  # package is not installed
    __version__ = "0.0.0"
