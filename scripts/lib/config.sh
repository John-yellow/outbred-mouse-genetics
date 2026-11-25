#!/usr/bin/env bash
set -euo pipefail
REPO_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
CONFIG_FILE=${PIPELINE_CONFIG:-"$REPO_ROOT/config/pipeline.env"}
if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "Missing config: $CONFIG_FILE (copy config/pipeline.env.example first)" >&2
  exit 2
fi
set -a
# shellcheck disable=SC1090
source "$CONFIG_FILE"
set +a
repo_path() {
  case "$1" in
    /*) printf '%s\n' "$1" ;;
    *) printf '%s/%s\n' "$REPO_ROOT" "$1" ;;
  esac
}
require_tool() {
  command -v "$1" >/dev/null 2>&1 || { echo "Required executable not found: $1" >&2; exit 2; }
}
