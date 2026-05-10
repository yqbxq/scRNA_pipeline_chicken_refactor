eda_gate_status() {
  local gate_id="$1"
  ensure_eda_control_files
  awk -F '\t' -v gate="${gate_id}" '
    NR > 1 && $1 == gate {
      print tolower($2)
      found = 1
      exit
    }
    END {
      if (!found) {
        print "pending"
      }
    }
  ' "${EDA_GATE_FILE}"
}

eda_gate_passed() {
  local gate_id="$1"
  local gate_status
  gate_status="$(eda_gate_status "${gate_id}")"
  [[ "${gate_status}" == "approved" || "${gate_status}" == "yes" || "${gate_status}" == "true" ]]
}

set_eda_gate_status() {
  local gate_id="$1"
  local gate_status="$2"
  local approved_by="${3:-}"
  local notes="${4:-}"

  ensure_eda_control_files

  local python_bin
  python_bin="$(detect_python)"

  env \
    EDA_GATE_FILE="${EDA_GATE_FILE}" \
    TARGET_GATE_ID="${gate_id}" \
    TARGET_GATE_STATUS="${gate_status}" \
    TARGET_APPROVED_BY="${approved_by}" \
    TARGET_NOTES="${notes}" \
    "${python_bin}" - <<'PY'
import csv
import os
from datetime import datetime
from pathlib import Path

gate_file = Path(os.environ["EDA_GATE_FILE"])
target_gate_id = os.environ["TARGET_GATE_ID"]
target_gate_status = os.environ["TARGET_GATE_STATUS"]
target_approved_by = os.environ["TARGET_APPROVED_BY"]
target_notes = os.environ["TARGET_NOTES"]

rows = []
with gate_file.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fieldnames = reader.fieldnames or ["gate_id", "status", "approved_by", "notes"]
    approved_field = "approved_by" if "approved_by" in fieldnames else ("reviewed_by" if "reviewed_by" in fieldnames else "approved_by")
    if approved_field not in fieldnames:
        fieldnames.append(approved_field)
    if "approved_at" not in fieldnames:
        fieldnames.append("approved_at")
    if "gate_profile" not in fieldnames:
        fieldnames.append("gate_profile")
    if "notes" not in fieldnames:
        fieldnames.append("notes")
    for row in reader:
        row.setdefault("approved_at", "")
        row.setdefault("gate_profile", "default")
        if row["gate_id"] == target_gate_id:
            row["status"] = target_gate_status
            row[approved_field] = target_approved_by
            if target_gate_status.lower() in {"approved", "yes", "true"} and not row.get("approved_at"):
                row["approved_at"] = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
            row["notes"] = target_notes
        rows.append(row)

if not any(row["gate_id"] == target_gate_id for row in rows):
    approved_at = ""
    if target_gate_status.lower() in {"approved", "yes", "true"}:
        approved_at = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    rows.append({
        "gate_id": target_gate_id,
        "status": target_gate_status,
        approved_field: target_approved_by,
        "approved_at": approved_at,
        "gate_profile": "default",
        "notes": target_notes,
    })

with gate_file.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    writer.writerows(rows)
PY
}
