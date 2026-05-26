#!/usr/bin/env python3
"""Build a project-level chicken-to-human ortholog LUT template.

This script is intentionally offline by default. It normalizes a user-supplied
TSV/CSV with chicken and human symbol columns into the communication LUT schema.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input TSV/CSV with chicken and human symbol columns.")
    parser.add_argument("--output", default="metadata/ortholog_chicken_human.tsv")
    parser.add_argument("--delimiter", default="\t")
    return parser.parse_args()


def pick(row: dict[str, str], names: list[str]) -> str:
    for name in names:
        value = row.get(name, "").strip()
        if value:
            return value
    return ""


def main() -> int:
    args = parse_args()
    in_path = Path(args.input)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with in_path.open("r", encoding="utf-8", newline="") as in_handle, out_path.open("w", encoding="utf-8", newline="") as out_handle:
        reader = csv.DictReader(in_handle, delimiter=args.delimiter)
        writer = csv.DictWriter(out_handle, delimiter="\t", fieldnames=["chicken_symbol", "human_symbol", "source", "confidence", "notes"])
        writer.writeheader()
        seen = set()
        for row in reader:
            chicken = pick(row, ["chicken_symbol", "chicken_gene", "source_symbol", "external_gene_name"])
            human = pick(row, ["human_symbol", "human_gene", "target_symbol", "hsapiens_homolog_associated_gene_name"])
            if not chicken or not human or (chicken.upper(), human.upper()) in seen:
                continue
            seen.add((chicken.upper(), human.upper()))
            writer.writerow({
                "chicken_symbol": chicken,
                "human_symbol": human,
                "source": pick(row, ["source", "orthology_type"]) or "user_supplied",
                "confidence": pick(row, ["confidence", "pair_quality"]) or "medium",
                "notes": pick(row, ["notes"]) or "",
            })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
