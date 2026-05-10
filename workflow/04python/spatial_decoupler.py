#!/usr/bin/env python3
from __future__ import annotations

import argparse


def main() -> int:
    parser = argparse.ArgumentParser(description="decoupleR-py spatial activity placeholder.")
    parser.add_argument("command", nargs="?", default="help")
    parser.parse_args()
    raise SystemExit("spatial_decoupler.py is a placeholder until ST decoupleR bridge logic is implemented.")


if __name__ == "__main__":
    raise SystemExit(main())
