#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT_DIR}"

Rscript - <<'RS'
source("workflow/05single_script/helpers/runtime_utils.R", encoding = "UTF-8")
source("workflow/05single_script/helpers/config.R", encoding = "UTF-8")
source("workflow/05single_script/helpers/ortholog_lookup_utils.R", encoding = "UTF-8")
source("workflow/05single_script/helpers/communication_consensus_utils.R", encoding = "UTF-8")
source("workflow/05single_script/helpers/multinichenet_wrapper.R", encoding = "UTF-8")

axis <- standardize_lr_axis_id("WNT5B_WNT5A", "FZD3+LRP6", "pGC", "rgGC")
stopifnot(identical(axis, "WNT5A|WNT5B|FZD3|LRP6|pGC->rgGC"))

df <- data.frame(
  ligand_human = "WNT5B_WNT5A",
  receptor_human = "FZD3+LRP6",
  source = "pGC",
  target = "rgGC",
  stringsAsFactors = FALSE
)
df <- standardize_lr_axis_id_in_df(df)
stopifnot(identical(df$lr_axis_id[[1]], axis))

available <- multinichenet_available_07c()
stopifnot(is.logical(available), length(available) == 1L)
RS

rg -q "bioconductor-multinichenetr" envs/environment_r_interaction.yml

echo "smoke_07c_multinichenet_ok"
