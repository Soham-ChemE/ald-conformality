"""Shared paths for the test suite. The package is imported from src/ (PYTHONPATH=src or pip install -e .)."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for sub in ("src", "coupling"):
    p = str(ROOT / sub)
    if p not in sys.path:
        sys.path.insert(0, p)
