#!/usr/bin/env python3
"""Run (and analyse) one or all PICSP3D feature examples.

    python3 run_example.py                       # list examples
    python3 run_example.py 02_two_stream_instability
    python3 run_example.py all                   # everything (takes a while)
    python3 run_example.py 05_dipole_trapping --analyze-only
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def examples():
    return sorted(d for d in os.listdir(HERE) if d[:2].isdigit() and os.path.isdir(os.path.join(HERE, d)))


def run(name, analyze_only=False):
    d = os.path.join(HERE, name)
    spec = importlib.util.spec_from_file_location(name, os.path.join(d, "analyze.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    if not analyze_only:
        mod.simulate_example(d)
    mod.analyze(d)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    only = "--analyze-only" in sys.argv
    if not args:
        print("\n".join(examples()))
        sys.exit(0)
    for n in examples() if args[0] == "all" else args:
        print(f"===== {n}")
        run(n.rstrip("/"), only)
