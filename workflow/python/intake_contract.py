from __future__ import annotations

import gzip
import re
from pathlib import Path

FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
ALLOWED_PLATFORMS = {"10x_cellranger", "dnbelab_c", "generic_mex", "auto", ""}
ALLOWED_GENE_ID_TYPES = {"symbol", "ensembl", "mixed", "unknown", "auto", ""}

ENSEMBL_GENE_RE = re.compile(r"^ENS[A-Z0-9]*G[0-9]+(?:\.[0-9]+)?$")
SYMBOL_RE = re.compile(r"^[A-Za-z][A-Za-z0-9._-]*$")


def normalize_value(raw: str | None) -> str:
    return (raw or "").strip()


def normalize_flag(raw: str | None, default: str = "yes") -> str:
    value = normalize_value(raw).lower()
    return value or default


def normalize_platform(raw: str | None) -> str:
    value = normalize_value(raw).lower()
    return value if value in ALLOWED_PLATFORMS else ""


def normalize_gene_id_type(raw: str | None) -> str:
    value = normalize_value(raw).lower()
    return value if value in ALLOWED_GENE_ID_TYPES else ""


def path_or_empty(path: Path | None) -> str:
    return str(path) if path is not None else ""


def has_fastq_files(path: Path) -> bool:
    if not path.exists() or not path.is_dir():
        return False
    for child in path.iterdir():
        if child.is_file() and child.name.endswith(FASTQ_SUFFIXES):
            return True
    return False


def is_mex_matrix_dir(path: Path) -> bool:
    if not path.exists() or not path.is_dir():
        return False
    matrix = (path / "matrix.mtx.gz").exists() or (path / "matrix.mtx").exists()
    features = any((path / name).exists() for name in ("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"))
    barcodes = (path / "barcodes.tsv.gz").exists() or (path / "barcodes.tsv").exists()
    return matrix and features and barcodes


def unique_paths(paths: list[Path]) -> list[Path]:
    ordered: list[Path] = []
    seen: set[str] = set()
    for path in paths:
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        ordered.append(path)
    return ordered


def first_existing_path(candidates: list[Path], expect_dir: bool = False) -> Path | None:
    for candidate in unique_paths(candidates):
        if expect_dir:
            if candidate.exists() and candidate.is_dir():
                return candidate
        elif candidate.exists():
            return candidate
    return None


def first_mex_dir(candidates: list[Path]) -> Path | None:
    for candidate in unique_paths(candidates):
        if is_mex_matrix_dir(candidate):
            return candidate
    return None


def resolve_cellranger_dir(sample_id: str, input_mode: str, source_path: Path, cellranger_out_dir: Path) -> Path | None:
    candidates: list[Path] = []
    if input_mode == "fastq":
        candidates.append(cellranger_out_dir / sample_id)
    elif input_mode == "cellranger_out":
        candidates.extend((source_path, source_path / sample_id, cellranger_out_dir / sample_id))
    elif input_mode == "matrix":
        candidates.extend((source_path, source_path / sample_id))
    for candidate in unique_paths(candidates):
        if (candidate / "outs").exists():
            return candidate
    return None


def resolve_dnbelab_root(sample_id: str, source_path: Path, data_dir: Path) -> Path | None:
    candidates = [
        source_path,
        source_path / sample_id,
        data_dir / sample_id,
        (data_dir / sample_id).parent / sample_id,
    ]
    for candidate in unique_paths(candidates):
        if (candidate / "filter_matrix").exists() or (candidate / "raw_matrix").exists():
            return candidate
    return None


def resolve_generic_mex_dir(sample_id: str, source_path: Path, data_dir: Path) -> Path | None:
    candidates = [
        source_path,
        source_path / "filtered_feature_bc_matrix",
        source_path / "filter_matrix",
        source_path / sample_id,
        source_path / sample_id / "filtered_feature_bc_matrix",
        source_path / sample_id / "filter_matrix",
        data_dir / sample_id,
    ]
    return first_mex_dir(candidates)


def detect_platform(sample_id: str, input_mode: str, source_path: Path, declared_platform: str, data_dir: Path, cellranger_out_dir: Path) -> str:
    if declared_platform in {"10x_cellranger", "dnbelab_c", "generic_mex"}:
        return declared_platform

    cellranger_dir = resolve_cellranger_dir(sample_id, input_mode, source_path, cellranger_out_dir)
    if cellranger_dir is not None:
        return "10x_cellranger"

    dnbelab_root = resolve_dnbelab_root(sample_id, source_path, data_dir)
    if dnbelab_root is not None:
        return "dnbelab_c"

    generic_dir = resolve_generic_mex_dir(sample_id, source_path, data_dir)
    if generic_dir is not None:
        return "generic_mex"

    if input_mode in {"fastq", "cellranger_out"}:
        return "10x_cellranger"
    if input_mode == "matrix":
        return "generic_mex"
    return "generic_mex"


def iter_features(feature_path: Path):
    opener = gzip.open if feature_path.suffix == ".gz" else open
    with opener(feature_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            yield line.split("\t")


def profile_feature_names(matrix_dir: Path | None, max_rows: int = 2000) -> dict[str, object]:
    empty_profile = {
        "feature_name_profile": "unknown",
        "feature_name_source": "",
        "n_features_checked": 0,
        "n_symbol_like": 0,
        "n_ensembl_like": 0,
    }
    if matrix_dir is None:
        return empty_profile

    feature_path = first_existing_path(
        [
            matrix_dir / "features.tsv.gz",
            matrix_dir / "features.tsv",
            matrix_dir / "genes.tsv.gz",
            matrix_dir / "genes.tsv",
        ]
    )
    if feature_path is None:
        return empty_profile

    total = 0
    symbol_like = 0
    ensembl_like = 0
    source = "matrix_features_col1"
    for fields in iter_features(feature_path):
        if len(fields) >= 2 and normalize_value(fields[1]):
            feature_name = fields[1].strip()
            source = "matrix_features_col2"
        else:
            feature_name = fields[0].strip()
            source = "matrix_features_col1"

        if not feature_name:
            continue

        total += 1
        if ENSEMBL_GENE_RE.match(feature_name):
            ensembl_like += 1
        elif SYMBOL_RE.match(feature_name):
            symbol_like += 1

        if total >= max_rows:
            break

    profile = "unknown"
    if total > 0:
        symbol_fraction = symbol_like / total
        ensembl_fraction = ensembl_like / total
        if symbol_fraction >= 0.80 and ensembl_fraction <= 0.10:
            profile = "symbol"
        elif ensembl_fraction >= 0.80 and symbol_fraction <= 0.10:
            profile = "ensembl"
        elif symbol_like > 0 and ensembl_like > 0:
            profile = "mixed"

    return {
        "feature_name_profile": profile,
        "feature_name_source": source,
        "n_features_checked": total,
        "n_symbol_like": symbol_like,
        "n_ensembl_like": ensembl_like,
    }


def resolve_gene_id_type(declared_gene_id_type: str, feature_profile: str) -> str:
    if declared_gene_id_type in {"symbol", "ensembl", "mixed", "unknown"}:
        return declared_gene_id_type
    if feature_profile in {"symbol", "ensembl", "mixed", "unknown"}:
        return feature_profile
    return "unknown"


def find_dnbelab_report(root: Path | None) -> Path | None:
    if root is None or not root.exists():
        return None
    reports = sorted(root.glob("*_scRNA_report.html"))
    return reports[0] if reports else None


def resolve_sample_contract(
    row: dict[str, str],
    fastq_dir: Path,
    cellranger_out_dir: Path,
    data_dir: Path,
) -> dict[str, object]:
    sample_id = normalize_value(row.get("sample_id"))
    input_mode = normalize_value(row.get("input_mode"))
    source_path = Path(normalize_value(row.get("source_path")) or ".")
    declared_platform = normalize_platform(row.get("platform"))
    declared_gene_id_type = normalize_gene_id_type(row.get("gene_id_type"))

    platform_resolved = detect_platform(
        sample_id=sample_id,
        input_mode=input_mode,
        source_path=source_path,
        declared_platform=declared_platform,
        data_dir=data_dir,
        cellranger_out_dir=cellranger_out_dir,
    )

    fastq_candidates: list[Path] = []
    if input_mode == "fastq":
        fastq_candidates.extend((source_path, source_path / sample_id, fastq_dir))
    has_fastq = any(has_fastq_files(candidate) for candidate in unique_paths(fastq_candidates))

    cellranger_dir = resolve_cellranger_dir(sample_id, input_mode, source_path, cellranger_out_dir)
    dnbelab_root = resolve_dnbelab_root(sample_id, source_path, data_dir)
    generic_filtered = resolve_generic_mex_dir(sample_id, source_path, data_dir)

    resolved_input_root: Path | None = None
    filtered_matrix_dir: Path | None = None
    raw_matrix_dir: Path | None = None
    metrics_path: Path | None = None
    bam_path: Path | None = None
    web_summary_path: Path | None = None
    raw_matrix_kind = "missing"

    if platform_resolved == "10x_cellranger":
        resolved_input_root = cellranger_dir or source_path
        if cellranger_dir is not None:
            outs_dir = cellranger_dir / "outs"
            filtered_matrix_dir = first_mex_dir([outs_dir / "filtered_feature_bc_matrix"])
            raw_matrix_dir = first_mex_dir([outs_dir / "raw_feature_bc_matrix"])
            raw_h5 = first_existing_path([outs_dir / "raw_feature_bc_matrix.h5"])
            if raw_matrix_dir is not None:
                raw_matrix_kind = "mex_dir"
            elif raw_h5 is not None:
                raw_matrix_dir = raw_h5
                raw_matrix_kind = "h5_file"
            metrics_path = first_existing_path([outs_dir / "metrics_summary.csv", cellranger_dir / "metrics_summary.csv"])
            bam_path = first_existing_path([outs_dir / "possorted_genome_bam.bam"])
            web_summary_path = first_existing_path([outs_dir / "web_summary.html"])
        if filtered_matrix_dir is None:
            filtered_matrix_dir = first_mex_dir(
                [
                    source_path,
                    source_path / "filtered_feature_bc_matrix",
                    source_path / sample_id,
                    data_dir / sample_id,
                ]
            )
        if filtered_matrix_dir is not None and resolved_input_root is None:
            resolved_input_root = filtered_matrix_dir.parent
    elif platform_resolved == "dnbelab_c":
        resolved_input_root = dnbelab_root or source_path
        if resolved_input_root is not None:
            filtered_matrix_dir = first_mex_dir([resolved_input_root / "filter_matrix"])
            raw_matrix_dir = first_mex_dir([resolved_input_root / "raw_matrix"])
            raw_matrix_kind = "mex_dir" if raw_matrix_dir is not None else "missing"
            metrics_path = first_existing_path(
                [
                    resolved_input_root / "metrics_summary.xls",
                    resolved_input_root / "metrics_summary.tsv",
                    resolved_input_root / "metrics_summary.csv",
                ]
            )
            bam_path = first_existing_path([resolved_input_root / "anno_decon_sorted.bam"])
            web_summary_path = find_dnbelab_report(resolved_input_root)
    else:
        filtered_matrix_dir = generic_filtered
        resolved_input_root = filtered_matrix_dir.parent if filtered_matrix_dir is not None else source_path
        raw_matrix_dir = first_mex_dir([source_path / "raw_matrix", source_path / "raw_feature_bc_matrix"])
        if raw_matrix_dir is not None:
            raw_matrix_kind = "mex_dir"
        metrics_path = first_existing_path(
            [
                source_path / "metrics_summary.tsv",
                source_path / "metrics_summary.csv",
                source_path / "metrics_summary.xls",
            ]
        )
        bam_path = first_existing_path(
            [
                source_path / "possorted_genome_bam.bam",
                source_path / "anno_decon_sorted.bam",
            ]
        )
        web_summary_path = first_existing_path([source_path / "web_summary.html"])

    profile_target = filtered_matrix_dir if filtered_matrix_dir is not None and filtered_matrix_dir.is_dir() else None
    if profile_target is None and raw_matrix_dir is not None and raw_matrix_dir.is_dir():
        profile_target = raw_matrix_dir
    feature_profile = profile_feature_names(profile_target)
    gene_id_type_resolved = resolve_gene_id_type(
        declared_gene_id_type=declared_gene_id_type,
        feature_profile=str(feature_profile["feature_name_profile"]),
    )

    has_filtered_matrix = filtered_matrix_dir is not None
    has_raw_matrix = raw_matrix_dir is not None
    has_bam = bam_path is not None
    has_metrics_summary = metrics_path is not None
    has_web_summary = web_summary_path is not None
    has_cellranger_out = cellranger_dir is not None

    return {
        "sample_id": sample_id,
        "input_mode": input_mode,
        "source_path": str(source_path),
        "platform": declared_platform or "auto",
        "platform_resolved": platform_resolved,
        "gene_id_type": declared_gene_id_type or "auto",
        "gene_id_type_resolved": gene_id_type_resolved,
        "reference_version": normalize_value(row.get("reference_version")),
        "resolved_input_root": path_or_empty(resolved_input_root),
        "cellranger_sample_dir": path_or_empty(cellranger_dir),
        "filtered_matrix_dir": path_or_empty(filtered_matrix_dir),
        "raw_matrix_dir": path_or_empty(raw_matrix_dir),
        "raw_matrix_kind": raw_matrix_kind,
        "metrics_path": path_or_empty(metrics_path),
        "bam_path": path_or_empty(bam_path),
        "web_summary_path": path_or_empty(web_summary_path),
        "has_fastq": str(has_fastq).lower(),
        "has_cellranger_out": str(has_cellranger_out).lower(),
        "has_filtered_matrix": str(has_filtered_matrix).lower(),
        "has_raw_matrix": str(has_raw_matrix).lower(),
        "has_bam": str(has_bam).lower(),
        "has_web_summary": str(has_web_summary).lower(),
        "has_metrics_summary": str(has_metrics_summary).lower(),
        "feature_name_profile": str(feature_profile["feature_name_profile"]),
        "feature_name_source": str(feature_profile["feature_name_source"]),
        "n_features_checked": int(feature_profile["n_features_checked"]),
        "n_symbol_like": int(feature_profile["n_symbol_like"]),
        "n_ensembl_like": int(feature_profile["n_ensembl_like"]),
    }
