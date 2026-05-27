# Communication Method Landscape

Module 07 now reserves one numbered slot per communication evidence layer:

| Stage | Method family | Role | Status |
|---|---|---|---|
| 07a | CellChat | Fast ligand-receptor hypothesis generation; never primary by itself. | implemented, `hypothesis_only` |
| 07b | LIANA+ consensus | Python-side multi-method ligand-receptor consensus. | implemented with dependency-aware empty-output fallback |
| 07c | NicheNet / MultiNicheNet | Receiver program-aware ligand prioritization. | legacy NicheNet executor with MultiNicheNet mode detection and schema hooks |
| 07d | Communication consensus | Cross-method evidence-tier assignment. | implemented |
| 07e | Communication EDA | Four-panel evidence report and review layer. | implemented |

07a writes `method_evidence_class=hypothesis_only` and `can_be_primary=no` into its task index and LR tables. CellChat-only axes are reportable candidate hypotheses, but final primary evidence requires support from LIANA, NicheNet/MultiNicheNet, or later COMMOT spatial evidence.

07b reads the scRNA H5AD mirror, filters communication labels through the 04c cluster eligibility table, and writes `liana_consensus_lr.tsv` with a canonical `lr_axis_id`. If LIANA/anndata are unavailable and `COMMUNICATION_REQUIRE_FULL_PIPELINE=no`, 07b writes an empty schema-valid table and manifest so later consensus stages can distinguish runtime absence from biological absence.

The current 07 chain still keeps `COMMUNICATION_REQUIRE_FULL_PIPELINE=no` until MultiNicheNet, CellChat downgrading, consensus, EDA, and COMMOT are all landed.

07c defaults to `NICHENET_MODE=auto`. When `multinichenetr` and `SingleCellExperiment` are unavailable, it records `multinichenet_mode_used=nichenet_legacy` and keeps the existing NicheNet executor. When the full MultiNicheNet runtime and fixtures are approved, the same mode field lets 07e/07d report whether evidence came from multi-sample MultiNicheNet or the single-condition fallback.

07d writes `method_consensus.tsv` and `evidence_tiers.tsv`. It assigns four reporting tiers: `primary`, `exploratory`, `candidate`, and `blocked`. The primary tier requires LIANA consensus and, by default, NicheNet/MultiNicheNet downstream support.

07e reads the 07d consensus table as the only authority for reportable axes. It
writes the 4-panel evidence report, a filtered axis table, panel status TSV, and
an optional HTML companion. Pair-specific report requirements come from
`communication_pairs.tsv` columns such as `methods_required`,
`min_methods_agreed`, `require_downstream_de`, `require_spatial_support`, and
`evidence_tier_required`.

Out-of-scope methods for this branch include CellNEST, DeepTalk, Tensor-cell2cell, and scTenifoldKnk. They can be added later only if they improve the evidence chain beyond the CellChat, LIANA+, MultiNicheNet, consensus, and COMMOT stack.
