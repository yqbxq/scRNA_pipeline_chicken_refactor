#!/usr/bin/env python3
"""修复 conda 的 bioconductor-data-packages 下载脚本。

这个脚本会扫描 conda 包缓存和项目环境中的 installBiocDataPackage.sh，
统一替换成一个更稳的版本：
1. 优先使用 depot.galaxyproject.org。
2. 给 curl 增加失败重试、连接超时和低速退出。
3. 支持从 BIOC_DATA_CACHE_DIR 读取预先准备好的 tar.gz，避免服务器重复慢速下载。
"""

from __future__ import annotations

import os
import sys
from pathlib import Path


SCRIPT_TEMPLATE = """#!/bin/bash
set -ex
# takes a single parameter, the package name

SCRIPT_DIR="$( dirname -- "${BASH_SOURCE[0]}" )/../share/bioconductor-data-packages"
json="${SCRIPT_DIR}/dataURLs.json"
FN=$(yq ".\\"$1\\".fn" "${json}" | tr -d '"')
while IFS= read -r value; do
  URLS+=($value);
done < <(yq ".\\"$1\\".urls[]" "${json}")
PREFERRED_URLS=()
FALLBACK_URLS=()
for URL in "${URLS[@]}"; do
  CLEAN_URL=$(echo $URL | tr -d '"')
  if [[ "$CLEAN_URL" == *"depot.galaxyproject.org"* ]]; then
    PREFERRED_URLS+=($URL)
  else
    FALLBACK_URLS+=($URL)
  fi
done
URLS=("${PREFERRED_URLS[@]}" "${FALLBACK_URLS[@]}")
MD5=$(yq ".\\"$1\\".md5" "${json}" | tr -d '"')

STAGING=$PREFIX/share/"$1"
mkdir -p $STAGING
TARBALL=$STAGING/$FN
SUCCESS=0

verify_md5() {
  local candidate="$1"
  if [[ $(uname -s) == "Linux" ]]; then
    md5sum -c <<<"$MD5  $candidate"
  else
    [[ $(md5 $candidate | cut -f4 -d " ") == "$MD5" ]]
  fi
}

CACHE_DIR="${BIOC_DATA_CACHE_DIR:-$HOME/.cache/bioconductor-data-packages}"
CACHE_TARBALL=$CACHE_DIR/$1/$FN
if [[ -f "$CACHE_TARBALL" ]]; then
  cp "$CACHE_TARBALL" "$TARBALL"
  if verify_md5 "$TARBALL"; then
    SUCCESS=1
  else
    rm -f "$TARBALL"
  fi
fi

if [[ $SUCCESS != 1 ]]; then
  for URL in ${URLS[@]}; do
    URL=$(echo $URL | tr -d '"')
    curl -L --fail --retry 3 --connect-timeout 15 --speed-time 30 --speed-limit 2048 $URL > $TARBALL
    [[ $? == 0 ]] || continue
    if verify_md5 "$TARBALL"; then
      SUCCESS=1
      mkdir -p "$(dirname "$CACHE_TARBALL")"
      cp "$TARBALL" "$CACHE_TARBALL"
      break
    fi
  done
fi

if [[ $SUCCESS != 1 ]]; then
  echo "ERROR: post-link.sh was unable to download any of the following URLs with the md5sum $MD5:"
  printf '%s\\n' "${URLS[@]}"
  exit 1
fi

R CMD INSTALL --library=$PREFIX/lib/R/library $TARBALL
rm $TARBALL
rmdir $STAGING
"""


def iter_targets() -> list[Path]:
    targets: list[Path] = []

    conda_root_env = os.environ.get("CONDA_ROOT", "")
    conda_root = Path(conda_root_env) if conda_root_env else None
    project_root_env = os.environ.get("PROJECT_ROOT", "")
    project_root = Path(project_root_env) if project_root_env else None

    search_roots = []
    if conda_root:
        search_roots.append(conda_root / "pkgs")
    if project_root:
        search_roots.append(project_root / "envs" / "conda")
    if len(sys.argv) > 1:
        search_roots.extend(Path(arg) for arg in sys.argv[1:])

    seen: set[Path] = set()
    for root in search_roots:
        if not root.exists():
            continue
        for path in root.rglob("installBiocDataPackage.sh"):
            if path not in seen:
                targets.append(path)
                seen.add(path)
    return targets


def patch_file(path: Path) -> bool:
    current = path.read_text(encoding="utf-8")
    if current == SCRIPT_TEMPLATE:
        return False
    path.write_text(SCRIPT_TEMPLATE, encoding="utf-8")
    return True


def main() -> int:
    targets = iter_targets()
    if not targets:
        print("没有找到 installBiocDataPackage.sh，可在第一次 conda 下载后再运行。")
        return 0

    patched = [str(path) for path in targets if patch_file(path)]
    if patched:
        print("\n".join(patched))
    else:
        print("未发现需要修改的 bioconductor-data-packages 脚本。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
