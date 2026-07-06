# Project-local single scripts

This directory stores one-off, project-local runners that were used to produce
auditable analysis artifacts outside the formal numbered workflow stages.

Current scripts:

- `pre_sub_manual_panorama_ovary_annotation.R`
  - Applies the user-curated `Panorama ovary annotation` / `卵巢全景注释结果`
    mapping to `/home/user_test/syf_f5/pre-sub/F6_1_panorama_project_panel_v3`.
  - Reads `results/checkpoints/04_after_subcluster_annotation.rds`.
  - Writes manual annotation tables, UMAPs, reports, and
    `results/checkpoints/05_after_manual_panorama_ovary_annotation.rds`.

- `pre_sub_gc_manual_subcluster.R`
  - Subclusters manually annotated GC-like / follicular somatic cells from the
    same panel v3 project.
  - Reads `results/checkpoints/05_after_manual_panorama_ovary_annotation.rds`.
  - Writes GC manual subcluster marker/evidence tables, UMAPs, reports, and
    `results/checkpoints/06_after_manual_GC_subcluster_annotation.rds`.

These scripts intentionally default to the server project paths used for the
F6_1 pre-sub analysis. Override `PROJECT_ROOT` and related environment
variables before reuse on a different project.
