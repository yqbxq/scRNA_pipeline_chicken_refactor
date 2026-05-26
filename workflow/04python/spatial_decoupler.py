#!/usr/bin/env python3
from __future__ import annotations

import argparse

from helpers.spatial_io import resolve_spatial_h5ad_path


def main() -> int:
    parser = argparse.ArgumentParser(description="decoupleR-py spatial activity placeholder.")
    parser.add_argument("command", nargs="?", default="help")
    parser.add_argument("--st-h5ad")
    parser.add_argument("--st-module", default="spatial_05_marker")
    args = parser.parse_args()
    if args.command == "resolve-h5ad":
        print(resolve_spatial_h5ad_path(module=args.st_module, fallback_path=args.st_h5ad))
        return 0
    raise SystemExit("spatial_decoupler.py is a placeholder until ST decoupleR bridge logic is implemented.")


if __name__ == "__main__":
    raise SystemExit(main())
