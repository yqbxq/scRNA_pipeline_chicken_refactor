#!/usr/bin/env bash

# 兼容入口：真实配置已经迁移到 config/server_config.sh。
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "${SCRIPT_DIR}/config/server_config.sh"
