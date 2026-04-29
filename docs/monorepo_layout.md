# Monorepo Layout

This branch is the integration branch for the single-repository pipeline layout.
The active branch name is `all-repo`.

## Directory Contract

```text
scRNA_pipeline_chicken_refactor/
├── config/                     Project and marker-panel configuration
├── docs/                       Human-readable design and migration notes
├── envs/                       Conda and post-install environment recipes
├── metadata/                   Project metadata and generated module inputs
├── workflow/                   Canonical pipeline implementation
│   ├── 01run.sh                Stage dispatcher
│   ├── 02lib/                  Shared Bash runtime, gates, stage helpers
│   ├── 03stages/               Stage wrappers, numbered by analysis module
│   ├── 04python/               Workflow-owned Python helpers
│   ├── 05single_script/        Workflow-owned single-module R/Bash scripts
│   ├── 06tools/                Small utility scripts
│   └── testing/                Smoke tests for active workflow code
```

## Active Development Rule

The full tree lives together on `all-repo`. Feature work should branch from this
branch and keep the same directory layout:

- `feat/m1-analysis-questions`
- `feat/m2-metadata-validator`
- `feat/m3-metadata-generator`
- `feat/m4-stage-metadata-hook`
- `feat/m5-cell-subtype-backfill`

Do not use folder-per-branch for active development. A Git branch represents a
whole repository tree, not one mounted subdirectory.

## Metadata Direction

The intended metadata architecture is:

- Tier 0: static experimental design tables, such as `metadata/samples.tsv`
- Tier 1: `metadata/analysis_questions.tsv`, the human-maintained source of truth
- Tier 2: generated module input tables, not hand-edited
- Tier 3: `workflow/03stages` wrappers and `workflow/05single_script` modules

The placeholder files in `metadata/` reserve the planned paths. M1-M5 should
replace placeholders with the implemented schema, validator, generator, stage
hook, and cell-subtype backfill.
