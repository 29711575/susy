#!/usr/bin/env python3
"""Compatibility entry; the unchanged static implementation now lives in limit/source."""
import sys
from pathlib import Path

__file__ = str(Path(__file__).resolve().parents[1] / "limit" / "source" / Path(__file__).name)
sys.path.insert(0, str(Path(__file__).parent))
exec(compile(Path(__file__).read_bytes(), __file__, "exec"), globals(), globals())
