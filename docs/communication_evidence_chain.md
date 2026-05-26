# Communication Evidence Chain

The 07 module is organized as a five-stage evidence chain:

| Evidence layer | Stage | Output intent |
|---|---|---|
| Hypothesis | 07a CellChat | Broad ligand-receptor candidates and pathway hypotheses. |
| LR consensus | 07b LIANA+ | Agreement across ligand-receptor methods. |
| Receiver program support | 07c NicheNet/MultiNicheNet | Ligands linked to receiver gene programs. |
| Evidence tier | 07d consensus | A single per-axis tier from cross-method support. |
| Review | 07e EDA | Tables, plots, and gate-facing summary. |

The planned evidence tiers are reserved but not fully assigned in P-R03-A:

| Tier | Intended meaning |
|---|---|
| `primary` | Supported by consensus LR evidence and receiver program evidence. |
| `exploratory` | Supported by one major method family or limited cell-count evidence. |
| `candidate_only` | Hypothesis-level communication requiring validation. |
| `skip` | Insufficient inventory, DEG, or method support for interpretation. |

P-R03-A only creates the stable stage names and placeholder manifests. The scoring and tier decision matrix are implemented in later P-R03 PRs.
