# Legacy Branch Migration Map

This branch consolidates the existing milestone history into the active
single-repository layout.

## Source Branches

| Source branch | Source head | Imported into active layout |
|---|---:|---|
| `main` | `5298f2d` | Original baseline, retained in repository history |
| `v00` | `fb34eea` | `shell/02lib`, `shell/03stages/00-02`, `shell/05single_script/00-02` |
| `v01` | `f810930` | Archived standalone 01 attempt; superseded by shell-owned 01 module in `v00` lineage |
| `v03` | `29016d0` | `shell/03stages/03_panorama.sh`, `shell/05single_script/03*`, panorama helpers |
| `v04` | `447bfa2` | `shell/03stages/04*`, `shell/05single_script/04*`, subcluster helpers |
| `v05` | `441abd1` | `shell/03stages/05_deg.sh`, `shell/05single_script/05*`, DEG helpers |
| `v06` | `07fcc60` | `shell/03stages/06_enrichment.sh`, `shell/05single_script/06*`, enrichment helpers |
| `v07` | `e0d5208` | `shell/03stages/07_communication.sh`, `shell/05single_script/07*`, interaction helpers |
| `v08` | local WIP | Imported as commit `feat(08): import regulation WIP into all-repo` |

## History Handling

The `all-repo` branch starts from `origin/v07`, so the `v00` through `v07`
milestone commits are part of normal ancestry. The local `v08` worktree changes
were imported as a dedicated commit because they were not committed on `v08`.

The independent `v01` worktree branch is not replayed as active code because the
later shell-owned 01 implementation superseded it. Keep `origin/v01` as an
archive reference when investigating that attempt.

## Module Ownership

| Module | Active paths |
|---|---|
| 00 ortholog | `shell/03stages/00_ortholog.sh`, `shell/05single_script/00*` |
| 01 raw/pre-QC | `shell/03stages/01_build_raw.sh`, `shell/05single_script/01*` |
| 02 QC | `shell/03stages/02_qc.sh`, `shell/05single_script/02*` |
| 03 panorama | `shell/03stages/03_panorama.sh`, `shell/05single_script/03*` |
| 04 subcluster | `shell/03stages/04*`, `shell/05single_script/04*` |
| 05 DEG | `shell/03stages/05_deg.sh`, `shell/05single_script/05*` |
| 06 enrichment | `shell/03stages/06_enrichment.sh`, `shell/05single_script/06*` |
| 07 communication | `shell/03stages/07_communication.sh`, `shell/05single_script/07*` |
| 08 regulation | `shell/03stages/08_regulation.sh`, `shell/05single_script/08*` |
| 09 trajectory | `shell/03stages/09_trajectory.sh` |
| 10 velocity | `shell/03stages/10_velocity.sh`, `shell/04python/scvelo_pipeline.py` |

## Next Milestones

Use the following feature branches from `all-repo`:

- `feat/m1-analysis-questions`
- `feat/m2-metadata-validator`
- `feat/m3-metadata-generator`
- `feat/m4-stage-metadata-hook`
- `feat/m5-cell-subtype-backfill`
