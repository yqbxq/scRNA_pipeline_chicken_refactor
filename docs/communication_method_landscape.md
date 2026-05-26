# Communication Method Landscape

Module 07 now reserves one numbered slot per communication evidence layer:

| Stage | Method family | Role | Status |
|---|---|---|---|
| 07a | CellChat | Fast hypothesis generation from ligand-receptor probability networks. | implemented |
| 07b | LIANA+ consensus | Python-side multi-method ligand-receptor consensus. | placeholder for P-R03-B |
| 07c | NicheNet / MultiNicheNet | Receiver program-aware ligand prioritization. | existing NicheNet code renamed to 07c |
| 07d | Communication consensus | Cross-method evidence-tier assignment. | placeholder for P-R03-E |
| 07e | Communication EDA | Reporting and review layer. | renamed from old 07c |

The current P-R03-A change is a numbering and orchestration foundation only. It does not change CellChat or NicheNet scoring logic. `COMMUNICATION_REQUIRE_FULL_PIPELINE=no` lets placeholder stages write placeholder manifests while the later LIANA+, MultiNicheNet, consensus, and COMMOT PRs land.

Out-of-scope methods for this branch include CellNEST, DeepTalk, Tensor-cell2cell, and scTenifoldKnk. They can be added later only if they improve the evidence chain beyond the CellChat, LIANA+, MultiNicheNet, consensus, and COMMOT stack.
