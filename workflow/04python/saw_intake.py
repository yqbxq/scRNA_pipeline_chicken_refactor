#!/usr/bin/env python3
from __future__ import annotations

import argparse


def main() -> int:
    parser = argparse.ArgumentParser(description="SAW/STOmics bundle conversion placeholder.")
    parser.add_argument("--input", required=True, help="SAW/STOmics output directory or GEF file.")
    parser.add_argument("--output", required=True, help="Standardized data/<sample>/outs directory.")
    parser.parse_args()
    raise SystemExit(
        "saw_intake.py is a placeholder. M2 only audits SAW files; real GEF conversion is deferred to M5."
    )


if __name__ == "__main__":
    raise SystemExit(main())
