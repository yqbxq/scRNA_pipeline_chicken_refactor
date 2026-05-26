# Evidence Tier Policy

P-R01-A defines six evidence tiers for downstream DE and communication modules.

| Tier | Trigger | Downstream action |
|---|---|---|
| `primary` | Sample count, per-sample cells, total UMI, and single-sample dominance all pass. | Formal pseudobulk DE or formal communication analysis. |
| `exploratory` | Replicates are present, but a soft quantitative rule such as per-sample cell count or total UMI is weak. | Run with warning and report as exploratory. |
| `module_score_only` | Both groups exist but sample/cell support is too weak for DE and no parent merge is available. | Score gene programs only; no formal DE/communication. |
| `merge_to_parent` | Cluster is too small and a parent cluster is available. | Merge into parent and mark the result as parent-derived. |
| `candidate_only` | A single sample dominates the cluster. | Report as a hypothesis only; no formal DE/communication. |
| `skip` | One comparison group has zero usable samples or support falls below the soft floor. | Do not produce a biological result. |

The policy treats zero-sample group comparisons as non-comparable. Single-sample dominance is kept separate from ordinary exploratory analysis because it is a sample-specific hypothesis, not a replicated finding.

The first producer is 04c subcluster EDA, which emits `cluster_eligibility_tsv` for both scRNA and spatial subcluster layers. Later PRs consume those manifests in 05b/05a and communication modules; unknown future tiers should be handled conservatively as exploratory or skipped by downstream callers.
