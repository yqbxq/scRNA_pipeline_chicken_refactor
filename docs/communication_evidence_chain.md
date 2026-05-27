# Communication Evidence Chain

The 07 module is organized as a five-stage evidence chain:

| Evidence layer | Stage | Output intent |
|---|---|---|
| Hypothesis | 07a CellChat | Broad ligand-receptor candidates and pathway hypotheses. |
| LR consensus | 07b LIANA+ | Agreement across ligand-receptor methods. |
| Receiver program support | 07c NicheNet/MultiNicheNet | Ligands linked to receiver gene programs. |
| Evidence tier | 07d consensus | A single per-axis tier from cross-method support. |
| Review | 07e EDA | Tables, plots, and gate-facing summary. |

The assigned evidence tiers are:

| Tier | Meaning |
|---|---|
| `primary` | LIANA consensus plus the configured method count and downstream NicheNet/MultiNicheNet support when required. |
| `exploratory` | At least one non-CellChat evidence layer supports the axis, or CellChat has an independent supporting layer. |
| `candidate` | CellChat-only or other hypothesis-level evidence that needs independent validation. |
| `blocked` | No supporting communication evidence is available. |

07e consumes 07d `method_consensus.tsv` as the single authority for reporting.
CellChat-only axes remain `candidate` and cannot set `can_be_primary=yes`.
COMMOT spatial support is reserved for the ST-side R03-G integration and is
reported as absent until a real `commot_spatial_hit` column is populated.
