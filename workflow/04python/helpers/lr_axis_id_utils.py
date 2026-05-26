from __future__ import annotations

import re
from typing import Any


def _is_missing(value: Any) -> bool:
    if value is None:
        return True
    text = str(value).strip()
    return text == "" or text.lower() in {"nan", "na", "none", "null", "-"}


def sort_complex(name: Any, sep: str = r"[_+/]") -> str:
    if _is_missing(name):
        return ""
    parts = [part.strip() for part in re.split(sep, str(name)) if part.strip()]
    return "|".join(sorted(parts))


def standardize_lr_axis_id(ligand: Any, receptor: Any, sender_cell_type: Any, receiver_cell_type: Any) -> str:
    ligand_complex = sort_complex(ligand)
    receptor_complex = sort_complex(receptor)
    sender = "" if _is_missing(sender_cell_type) else str(sender_cell_type).strip()
    receiver = "" if _is_missing(receiver_cell_type) else str(receiver_cell_type).strip()
    return f"{ligand_complex}|{receptor_complex}|{sender}->{receiver}"


def standardize_lr_axis_id_df(
    df,
    ligand_col: str = "ligand_human",
    receptor_col: str = "receptor_human",
    sender_col: str = "source",
    receiver_col: str = "target",
    out_col: str = "lr_axis_id",
):
    out = df.copy()
    out[out_col] = [
        standardize_lr_axis_id(ligand, receptor, sender, receiver)
        for ligand, receptor, sender, receiver in zip(
            out[ligand_col],
            out[receptor_col],
            out[sender_col],
            out[receiver_col],
        )
    ]
    return out
