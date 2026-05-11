# 08 Regulation Workflow

08 regulation is the regulatory evidence layer. It complements 05 DEG, 06 enrichment, and 07 communication; it does not replace them.

## Modules

| Module | Purpose | Main runtime |
| --- | --- | --- |
| `08a_scenic_export.R` | Export chicken expression as human-symbol matrix for SCENIC. | `r_main` |
| `08b_scenic_grn.sh` | Run pySCENIC GRNBoost2. | `pyscenic` |
| `08c_scenic_regulons.R` | Build motif-backed regulons and AUCell matrix. | `r_scenic` |
| `08d_scenic_downstream.R` | Generate RSS, CSI, figures, and SCENIC index tables. | `r_scenic` |
| `08e_decoupler.R` | Run decoupleR TF and pathway activity from mapped human prior networks. | `r_decoupler` |
| `08f_regulation_eda.R` | Build the full or partial regulation report. | `r_main` |

## Species Mapping

SCENIC and decoupleR intentionally use different mapping directions:

| Method | Mapping strategy |
| --- | --- |
| SCENIC | Chicken expression matrix is mapped to human symbols before pySCENIC and cisTarget. |
| decoupleR | Human DoRothEA/PROGENy prior networks are mapped to chicken symbols; expression remains chicken gene space. |

This preserves chicken expression context for decoupleR while still using human curated priors.

## Running

Recommended first pass:

```bash
REGULATION_LAYERS=panorama bash workflow/03stages/08_regulation.sh
```

To add a subcluster:

```bash
REGULATION_LAYERS=panorama,GC_subcluster bash workflow/03stages/08_regulation.sh
```

decoupleR-only smoke:

```bash
REGULATION_LAYERS=panorama RUN_SCENIC_GRN=no RUN_DECOUPLER=yes bash workflow/03stages/08_regulation.sh
```

SCENIC resource-only preparation:

```bash
SCENIC_RESOURCES_ONLY=yes bash workflow/03stages/08_regulation.sh
```

## Partial Reports

`08f` is partial-report friendly:

| SCENIC | decoupleR | 08f behavior |
| --- | --- | --- |
| ok | ok | Full report. |
| ok | skipped | SCENIC-only partial report. |
| failed/skipped | ok | decoupleR-only partial report. |
| failed/skipped | failed/skipped | Stage fails unless `ALLOW_REGULATION_EMPTY_REPORT=yes`. |

05/06/07 evidence files are soft optional. Missing summaries are marked `not_available` in `module_status.tsv` and in the report.

## Interpretation

- SCENIC and decoupleR are complementary evidence streams, not equivalent methods.
- Overlap between SCENIC and decoupleR TFs strengthens confidence, but lack of overlap is not automatically an error.
- SCENIC depends on motif/regulon inference and human cisTarget resources.
- decoupleR depends on curated human prior networks mapped to chicken genes.
- 08 TF/pathway activity cannot replace DEG, enrichment, or communication evidence.
