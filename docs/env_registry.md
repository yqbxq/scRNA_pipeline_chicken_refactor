# Environment Registry

This file records project-level environment variables that are exported through `workflow/02lib/common.sh` and passed into runtime environments through `workflow/02lib/env_registry.sh`.

| Variable | Default | Purpose |
|---|---:|---|
| `CHECKPOINT_MODE` | `mtime` | Selects checkpoint behavior: `mtime`, `fingerprint`, or `off`. |
| `CHECKPOINT_FINGERPRINT_FALLBACK_ON_ERROR` | `mtime` | Controls fingerprint failure handling: warn and continue with mtime semantics, or `error`. |
| `CHECKPOINT_FINGERPRINT_VERBOSE` | `no` | Prints current and saved fingerprint values during stale checks when set to `yes`. |
| `CHECKPOINT_LARGE_FILE_THRESHOLD_MB` | `100` | Uses large-file mixed digests above this size for heavy binary inputs. |
| `MODULE_CHECKPOINT_FINGERPRINT_VERSION` | `1.0` | Version marker for the checkpoint fingerprint infrastructure. |
