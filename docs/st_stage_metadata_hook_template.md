# ST Stage Metadata Hook Template

Future ST and joint-analysis stages must consume generated Tier 2 metadata in
the same way as active scRNA stages.

After sourcing `workflow/02lib/common.sh`, call:

```bash
ensure_metadata_fresh
```

Recommended placement:

```bash
source "${WORKFLOW_ROOT}/02lib/common.sh"

ensure_metadata_fresh
check_stage_deps "11_st_qc"
```

`ensure_metadata_fresh` checks `metadata/analysis_questions.tsv` against all
generated Tier 2 metadata tables. If any table is missing or older than
`analysis_questions.tsv`, it runs `workflow/03stages/95_run_metadata_generator.sh`
and then validates with `workflow/03stages/96_validate_metadata.sh`.

Do not hand-edit generated Tier 2 files. Update
`metadata/analysis_questions.tsv`, then let the stage hook regenerate:

- `metadata/comparisons.tsv`
- `metadata/communication_pairs.tsv`
- `metadata/trajectory_pairs.tsv`
- `metadata/scenic_targets.tsv`
- `metadata/enrichment_targets.tsv`
- `metadata/deconv_pairs.tsv`
- `metadata/spatial_pairs.tsv`

Current ST rows in `analysis_questions.tsv` remain `status=planned`; ST-H owns
activating those rows and implementing the ST fan-out rules.
