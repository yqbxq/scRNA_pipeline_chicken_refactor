# Checkpoint Manifest Fingerprints

P-R02-A adds an optional `fingerprints` block to module `_manifest.json` files. The default checkpoint behavior remains `CHECKPOINT_MODE=mtime`, so existing runs keep the old mtime-only behavior until fingerprint mode is explicitly enabled.

## Schema

```json
{
  "fingerprints": {
    "algorithm": "sha256",
    "computed_at": "2026-05-25T12:00:00Z",
    "input_hash": "sha256-of-input-path-digests",
    "script_hash": "sha256-of-script-and-depth-1-helper-digests",
    "params_hash": "sha256-of-selected-env-params",
    "params_names": ["RANDOM_SEED"],
    "checkpoint_mode": "fingerprint"
  }
}
```

`input_hash` hashes each declared input path. Directories are expanded recursively with hidden files skipped. JSON manifests are expanded through their `outputs[*].path` entries so timestamps in upstream manifests do not create unstable downstream hashes.

`script_hash` hashes the stage script. R `source()` / `source_utf8()` and Python `helpers.*` imports are included at depth 1.

`params_hash` hashes only explicitly declared parameter names. Do not include secrets, tokens, or credentials in `params_names`.

## Decision Tree

`CHECKPOINT_MODE=mtime` uses the legacy `is_stale_output()` logic.

`CHECKPOINT_MODE=off` always reruns the stage.

`CHECKPOINT_MODE=fingerprint` compares current hashes to the saved manifest hashes:

1. If the manifest is missing or has no valid `fingerprints` block, the stage is stale.
2. If any of `input_hash`, `script_hash`, or `params_hash` differs, the stage is stale.
3. Otherwise the stage is skipped.
4. After a fingerprint-mode rerun, the stage runner writes a fresh `fingerprints` block back to the manifest.

If fingerprint writing fails and `CHECKPOINT_FINGERPRINT_FALLBACK_ON_ERROR=mtime`, the run continues with a warning. Set it to `error` to fail the stage.

## Upgrade Path

P-R02-A only installs the infrastructure. It does not change business stage call sites and does not switch the default mode. The first manual run with `CHECKPOINT_MODE=fingerprint` will mark old manifests stale once, write fingerprints, and then skip stable stages on later runs.
