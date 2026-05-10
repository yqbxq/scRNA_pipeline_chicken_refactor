#!/usr/bin/env python3
from __future__ import annotations

import argparse


def main() -> int:
    parser = argparse.ArgumentParser(description="cell2location spatial deconvolution placeholder.")
    parser.add_argument("command", nargs="?", default="help")
    parser.parse_args()
    raise SystemExit("cell2location_pipeline.py is a placeholder until cell2location deconvolution is implemented.")


if __name__ == "__main__":
    raise SystemExit(main())
