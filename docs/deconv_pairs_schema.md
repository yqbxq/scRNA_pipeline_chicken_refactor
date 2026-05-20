# Deconvolution Pairs Schema

`metadata/deconv_pairs.tsv` is generated from `metadata/analysis_questions.tsv` when ST fan-out is enabled.

| Column | Meaning |
| --- | --- |
| `deconv_id` | Stable deconvolution task identifier. |
| `source_question_id` | Source analysis question, usually `I08_deconv_panorama`. |
| `st_scope` | ST scope, currently `panorama_st`. |
| `reference_scope` | Reference subset; `panorama` or `all` are the default scopes. |
| `section_filter` | Section selector, with `*` or blank meaning all sections. |
| `condition_split_var` | Optional metadata column for condition-aware summaries. |
| `condition_split_values` | Optional comma-separated condition values. |
| `tool` | `all`, or a comma-separated subset of `rctd,transfer,card,cell2location`. |
| `enabled` | `yes` or `no`. |
| `notes` | Free-text provenance and review notes. |

Example:

```tsv
deconv_id	source_question_id	st_scope	reference_scope	section_filter	condition_split_var	condition_split_values	tool	enabled	notes
I08_deconv_panorama__panorama	I08_deconv_panorama	panorama_st	panorama	*	condition	syf,f5	all	yes	panorama reference deconvolution into ST sections
```

The current repository copy keeps `deconv_pairs.tsv` header-only until ST fan-out is enabled for a project.
