from __future__ import annotations

import hashlib
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable

LARGE_SUFFIXES = {".bam", ".cram", ".fastq", ".fq", ".gz", ".h5", ".h5ad", ".loom"}
MAX_MANIFEST_DEPTH = 5


def _sha256_bytes(chunks: Iterable[bytes]) -> str:
    digest = hashlib.sha256()
    for chunk in chunks:
        digest.update(chunk)
    return digest.hexdigest()


def compute_string_sha256(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def compute_file_sha256(path: str | Path, large_file_threshold_mb: int | None = None) -> str:
    item = Path(path)
    if not item.exists():
        return "missing"
    if item.is_dir():
        return "dir"

    threshold_mb = large_file_threshold_mb
    if threshold_mb is None:
        threshold_mb = int(os.environ.get("CHECKPOINT_LARGE_FILE_THRESHOLD_MB", "100") or "100")
    threshold_bytes = threshold_mb * 1024 * 1024
    stat = item.stat()
    if stat.st_size > threshold_bytes and {suffix.lower() for suffix in item.suffixes}.intersection(LARGE_SUFFIXES):
        with item.open("rb") as handle:
            head_hash = hashlib.sha256(handle.read(1024 * 1024)).hexdigest()
        return f"large:{stat.st_size}:{stat.st_mtime_ns}:{head_hash}"

    def chunks():
        with item.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                yield chunk

    return _sha256_bytes(chunks())


def _is_manifest(path: Path) -> bool:
    if not path.is_file() or path.suffix.lower() != ".json":
        return False
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return False
    return isinstance(data, dict) and isinstance(data.get("outputs"), dict)


def _expand_manifest_outputs(path: Path, depth: int, visited: set[Path]) -> list[Path]:
    if depth > MAX_MANIFEST_DEPTH:
        raise RuntimeError(f"manifest nesting exceeds {MAX_MANIFEST_DEPTH}: {path}")
    resolved = path.resolve()
    if resolved in visited:
        return []
    visited.add(resolved)

    data = json.loads(path.read_text(encoding="utf-8"))
    base_dir = Path(data.get("base_dir") or path.parent)
    output_paths: list[Path] = []
    for entry in (data.get("outputs") or {}).values():
        if not isinstance(entry, dict) or not entry.get("path"):
            continue
        candidate = Path(str(entry["path"]))
        if not candidate.is_absolute():
            candidate = base_dir / candidate
        output_paths.extend(_expand_path(candidate, depth + 1, visited))
    return output_paths


def _expand_path(path: Path, depth: int = 0, visited: set[Path] | None = None) -> list[Path]:
    if visited is None:
        visited = set()
    if not path.exists():
        return [path]
    if path.is_dir():
        files = []
        for child in path.rglob("*"):
            if not child.is_file():
                continue
            try:
                rel_parts = child.relative_to(path).parts
            except ValueError:
                rel_parts = child.parts
            if any(part.startswith(".") for part in rel_parts):
                continue
            files.append(child)
        return files
    if _is_manifest(path):
        return _expand_manifest_outputs(path, depth, visited)
    return [path]


def compute_input_hash(input_paths: Iterable[str | Path]) -> str:
    rows: list[tuple[str, str]] = []
    visited: set[Path] = set()
    for raw in input_paths:
        if str(raw) == "":
            continue
        for path in _expand_path(Path(raw), visited=visited):
            rows.append((str(path), compute_file_sha256(path)))
    aggregate = "".join(f"{path}\t{digest}\n" for path, digest in sorted(rows))
    return compute_string_sha256(aggregate)


def compute_script_hash(script_path: str | Path) -> str:
    script = Path(script_path)
    files = [script]
    if script.exists():
        content = script.read_text(encoding="utf-8", errors="ignore")
        if script.suffix == ".R":
            for match in re.finditer(r'(?:source|source_utf8)\s*\(\s*["\']([^"\']+)["\']\s*\)', content):
                helper = (script.parent / match.group(1)).resolve()
                if helper.exists() and helper.is_file():
                    files.append(helper)
        elif script.suffix == ".py":
            for match in re.finditer(r"(?:from\s+helpers\.(\w+)\s+import|import\s+helpers\.(\w+))", content):
                mod = match.group(1) or match.group(2)
                helper = (script.parent / "helpers" / f"{mod}.py").resolve()
                if helper.exists() and helper.is_file():
                    files.append(helper)

    rows = []
    seen: set[str] = set()
    for path in sorted(files, key=lambda item: str(item)):
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        rows.append((key, compute_file_sha256(path)))
    aggregate = "".join(f"{path}\t{digest}\n" for path, digest in rows)
    return compute_string_sha256(aggregate)


def compute_params_hash(params_names: Iterable[str], getenv: Callable[[str], str | None] | None = None) -> str:
    names = sorted({name.strip() for name in params_names if name and name.strip()})
    if not names:
        return ""
    if getenv is None:
        getenv = os.environ.get
    aggregate = "".join(f"{name}={getenv(name) or '__UNSET__'}\n" for name in names)
    return compute_string_sha256(aggregate)


def build_fingerprints(
    input_paths: Iterable[str | Path],
    script_path: str | Path,
    params_names: Iterable[str] = (),
    checkpoint_mode: str = "fingerprint",
) -> dict[str, object]:
    names = sorted({name.strip() for name in params_names if name and name.strip()})
    return {
        "algorithm": "sha256",
        "computed_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "input_hash": compute_input_hash(input_paths),
        "script_hash": compute_script_hash(script_path),
        "params_hash": compute_params_hash(names),
        "params_names": names,
        "checkpoint_mode": checkpoint_mode,
    }


def manifest_has_fingerprints(manifest_path: str | Path) -> bool:
    path = Path(manifest_path)
    if not path.exists():
        return False
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return False
    fp = data.get("fingerprints")
    return (
        isinstance(fp, dict)
        and fp.get("algorithm") == "sha256"
        and "input_hash" in fp
        and "script_hash" in fp
        and "params_hash" in fp
    )

