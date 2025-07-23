"""Lightweight function and line profiler integrated with faststat.

This module provides a :class:`Profiler` that attaches using the standard
``sys.setprofile`` and ``sys.settrace`` hooks.  A :class:`faststat.Duration`
instance is maintained for each function or line that executes, allowing
collection of high level timing statistics.
"""

from __future__ import annotations

import argparse
import runpy
import sys
from collections import defaultdict
from types import FrameType
from typing import Callable, Dict, Tuple

from .faststat import Duration, nanotime


class Profiler:
    """Collect execution timing statistics using faststat."""

    def __init__(self, *, line_stats: bool = False, func_stats: bool = True):
        self.line_stats = line_stats
        self.func_stats = func_stats
        self._call_start: Dict[FrameType, int] = {}
        self._line_info: Dict[FrameType, Tuple[int, int]] = {}
        self.func_data: Dict[Tuple[str, str, int], Duration] = defaultdict(
            lambda: Duration(interval=False)
        )
        self.line_data: Dict[Tuple[str, int], Duration] = defaultdict(
            lambda: Duration(interval=False)
        )

    # tracing ---------------------------------------------------------------
    def _trace(self, frame: FrameType, event: str, arg):
        now = nanotime()
        if event == "call":
            if self.func_stats:
                self._call_start[frame] = now
            if self.line_stats:
                self._line_info[frame] = (now, frame.f_lineno)
            return self._trace
        if event == "return":
            if self.func_stats and frame in self._call_start:
                start = self._call_start.pop(frame)
                key = (
                    frame.f_code.co_filename,
                    frame.f_code.co_name,
                    frame.f_code.co_firstlineno,
                )
                self.func_data[key].end(start)
            if self.line_stats and frame in self._line_info:
                start, lineno = self._line_info.pop(frame)
                self.line_data[(frame.f_code.co_filename, lineno)].end(start)
            return self._trace
        if event == "line" and self.line_stats:
            if frame in self._line_info:
                start, lineno = self._line_info[frame]
                self.line_data[(frame.f_code.co_filename, lineno)].end(start)
            self._line_info[frame] = (now, frame.f_lineno)
            return self._trace
        return self._trace

    # public API -----------------------------------------------------------
    def start(self) -> None:
        sys.setprofile(self._trace)
        if self.line_stats:
            sys.settrace(self._trace)

    def stop(self) -> None:
        sys.setprofile(None)
        if self.line_stats:
            sys.settrace(None)
        self._call_start.clear()
        self._line_info.clear()

    def report(self, *, limit: int = 20) -> None:
        """Print collected statistics sorted by mean duration."""

        if self.func_stats:
            print("Function timings:")
            for (file, name, lineno), stat in sorted(
                self.func_data.items(), key=lambda kv: kv[1].mean, reverse=True
            )[:limit]:
                print(
                    f"{name} {file}:{lineno} n={stat.n} mean={stat.mean} max={stat.max}"
                )
        if self.line_stats:
            print("Line timings:")
            for (file, lineno), stat in sorted(
                self.line_data.items(), key=lambda kv: kv[1].mean, reverse=True
            )[:limit]:
                print(
                    f"{file}:{lineno} n={stat.n} mean={stat.mean} max={stat.max}"
                )


# CLI ---------------------------------------------------------------------

def cli(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Profile python code using faststat")
    parser.add_argument(
        "-m",
        "--mode",
        choices=["function", "line", "both"],
        default="function",
        help="What to profile",
    )
    parser.add_argument("script", help="Python script to run")
    parser.add_argument("args", nargs=argparse.REMAINDER, help="Arguments for script")
    ns = parser.parse_args(argv)

    line_stats = ns.mode in ("line", "both")
    func_stats = ns.mode in ("function", "both")
    prof = Profiler(line_stats=line_stats, func_stats=func_stats)
    sys.argv = [ns.script] + ns.args
    prof.start()
    try:
        runpy.run_path(ns.script, run_name="__main__")
    finally:
        prof.stop()
    prof.report()


if __name__ == "__main__":  # pragma: no cover - manual usage
    cli()
