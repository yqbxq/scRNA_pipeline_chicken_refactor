from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any


class OrthologLookup:
    def __init__(self, path: str | Path | None = None):
        self.path = Path(path) if path else None
        self._lookup: dict[str, str] = {}
        if self.path and self.path.exists():
            self._load(self.path)

    def _load(self, path: Path) -> None:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                chicken = (
                    row.get("chicken_symbol")
                    or row.get("chicken_gene")
                    or row.get("source_symbol")
                    or row.get("external_gene_name")
                    or ""
                ).strip()
                human = (
                    row.get("human_symbol")
                    or row.get("human_gene")
                    or row.get("target_symbol")
                    or row.get("hsapiens_homolog_associated_gene_name")
                    or ""
                ).strip()
                if chicken and human and human != "-":
                    self._lookup[chicken.upper()] = human

    def chicken_to_human(self, symbol: Any) -> str:
        text = "" if symbol is None else str(symbol).strip()
        if not text or text == "-":
            return ""
        return self._lookup.get(text.upper(), text)

    def chicken_to_human_complex(self, value: Any) -> str:
        text = "" if value is None else str(value).strip()
        if not text or text == "-":
            return ""
        parts = [part.strip() for part in re.split(r"[_+/]", text) if part.strip()]
        mapped = [self.chicken_to_human(part) for part in parts]
        mapped = [part for part in mapped if part]
        return "_".join(mapped)
