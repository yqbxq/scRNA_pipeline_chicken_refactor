# Spatial Communication Pipeline

ST stage 08 verifies scRNA communication candidates in spatial context. COMMOT is one support layer, not a standalone mechanism call. Final spatial communication tiers combine scRNA consensus, ST deconvolution/validation, sender-receiver colocalization, 06d neighborhood support, and COMMOT signal.

## Stage Layout

| Stage | Script | Output |
|---|---|---|
| 08a | `workflow/05single_script/spatial/08a_spatial_communication_io.R` | `spatial_comm_candidates.tsv`, `spatial_comm_input_manifest.tsv`, `commot_input_manifest.tsv`, and `commot_lr_candidates.tsv`. |
| 08b Python | `workflow/04python/07f_commot_spatial.py` | `commot_lr.tsv`, one row per section and LR candidate. |
| 08b EDA | `workflow/05single_script/spatial/08b_spatial_communication_eda.R` | `commot_spatial_summary.tsv`, `commot_celltype_scores.tsv`, `commot_spot_pair_scores.tsv`, and a review report. |
| 08c | `workflow/05single_script/spatial/08c_lr_colocalization.R` | `lr_colocalization.tsv` and `sender_receiver_colocalization.tsv`. |
| 08d | `workflow/05single_script/spatial/08d_communication_neighborhood_consistency.R` | `communication_neighborhood_consistency.tsv`. |
| 08e | `workflow/05single_script/spatial/08e_spatial_communication_consensus.R` | `spatial_communication_consensus.tsv` and `spatial_communication_evidence_tier.tsv`. |
| 08f | `workflow/05single_script/spatial/08f_spatial_communication_report.R` | `spatial_communication_report.md`, `spatial_communication_question_gate_status.tsv`, and final panel candidates. |

`60_run_spatial_extensions.sh` runs this block after the deconvolution gate and
before niche derivation when `COMMOT_ENABLED != no`.

## Runtime Behavior

08b uses Python COMMOT when `commot` and `anndata` are available. If the runtime
is unavailable and `COMMOT_REQUIRE_RUNTIME=no`, the stage writes a schema-valid
skip table so downstream manifests remain explicit. Set
`COMMOT_REQUIRE_RUNTIME=yes` on the server when a missing COMMOT runtime should
block the run.

The COMMOT API used for the real runtime is
`commot.tl.cluster_communication_spatial_permutation`, which accepts an AnnData
object, an LR table, a distance threshold, and a cluster label column.

## 07d Integration

08b writes `results/spatial/tables/08_commot/commot_spatial_summary.tsv`.
07d reads this path through `COMMOT_SPATIAL_SUMMARY_TSV` and fills
`commot_spatial_hit` before assigning evidence tiers. 08e also reads the same
summary and crosses it with deconvolution, colocalization, and neighborhood
support. CellChat-only candidates remain `spatial_hypothesis`; `spatial_primary`
requires primary scRNA support, formal deconvolution support, sender/receiver
localization, neighborhood support, and either COMMOT support or ligand/receptor
expression-level colocalization. Sender/receiver abundance colocalization alone
is retained as localization evidence but does not replace COMMOT or LR expression
support for primary evidence.

Deconvolution support is tiered. `ok`/`PASS` can support primary evidence,
`warn`/`WARN` can support exploratory evidence, `ok_smoke` is smoke-only and
cannot produce `spatial_primary`, and missing/failed/skipped method outputs are
blocked.

## Question Gates

08f writes `spatial_communication_question_gate_status.tsv` for:

| Question | Meaning |
|---|---|
| `I19_comm_in_space` | Candidate communication has spatial support. |
| `I20_comm_neighbor_check` | Sender/receiver communication is consistent with 06d neighborhood evidence. |
| `I21_comm_final_figure` | Candidates are eligible for final spatial communication figure review. |
