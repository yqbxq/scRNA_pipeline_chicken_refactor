#!/usr/bin/env python3
from __future__ import annotations

import argparse


def main() -> int:
    parser = argparse.ArgumentParser(description="py_spatial execution placeholder for Squidpy/Scanpy ST tasks.")
    parser.add_argument("command", nargs="?", default="help")
    parser.parse_args()
    raise SystemExit("spatial_pipeline.py is a placeholder until the corresponding ST Python bridge stage is implemented.")


if __name__ == "__main__":
    raise SystemExit(main())
