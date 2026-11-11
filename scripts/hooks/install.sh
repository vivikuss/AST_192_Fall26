#!/usr/bin/env bash
# Install the blocking privacy hooks for this repository (DECISIONS E49):
#     scripts/hooks/install.sh
# Points core.hooksPath at scripts/hooks (repository-local config only).
set -euo pipefail
ROOT="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
chmod +x "$ROOT/scripts/hooks/pre-commit" "$ROOT/scripts/hooks/commit-msg" "$ROOT/scripts/hooks/pre-push"
git -C "$ROOT" config core.hooksPath scripts/hooks
echo "core.hooksPath = $(git -C "$ROOT" config core.hooksPath)"
"${AST192_HOOK_PYTHON:-python3}" -I "$ROOT/scripts/hooks/privacy_scan.py" staged
