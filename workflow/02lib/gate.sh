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
    for row in reader:
        if row["gate_id"] == target_gate_id:
            row["status"] = target_gate_status
            row["approved_by"] = target_approved_by
            row["notes"] = target_notes
        rows.append(row)

if not any(row["gate_id"] == target_gate_id for row in rows):
    rows.append({
        "gate_id": target_gate_id,
        "status": target_gate_status,
        "approved_by": target_approved_by,
        "notes": target_notes,
    })

with gate_file.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
PY
}
