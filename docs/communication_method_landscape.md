# Communication Method Landscape

Module 07 now reserves one numbered slot per communication evidence layer:

| Stage | Method family | Role | Status |
|---|---|---|---|
| 07a | CellChat | Fast hypothesis generation from ligand-receptor probability networks. | implemented |
| 07b | LIANA+ consensus | Python-side multi-method ligand-receptor consensus. | implemented with dependency-aware empty-output fallback |
| 07c | NicheNet / MultiNicheNet | Receiver program-aware ligand prioritization. | existing NicheNet code renamed to 07c |
| 07d | Communication consensus | Cross-method evidence-tier assignment. | placeholder for P-R03-E |
| 07e | Communication EDA | Reporting and review layer. | renamed from old 07c |

07b reads the scRNA H5AD mirror, filters communication labels through the 04c cluster eligibility table, and writes `liana_consensus_lr.tsv` with a canonical `lr_axis_id`. If LIANA/anndata are unavailable and `COMMUNICATION_REQUIRE_FULL_PIPELINE=no`, 07b writes an empty schema-valid table and manifest so later consensus stages can distinguish runtime absence from biological absence.

The current 07 chain still keeps `COMMUNICATION_REQUIRE_FULL_PIPELINE=no` until MultiNicheNet, CellChat downgrading, consensus, EDA, and COMMOT are all landed.

Out-of-scope methods for this branch include CellNEST, DeepTalk, Tensor-cell2cell, and scTenifoldKnk. They can be added later only if they improve the evidence chain beyond the CellChat, LIANA+, MultiNicheNet, consensus, and COMMOT stack.
