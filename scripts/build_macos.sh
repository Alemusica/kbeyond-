#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${MAX_SDK_ROOT:-}" ]]; then
    echo "ERROR: MAX_SDK_ROOT environment variable must point to the Max SDK" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="${SCRIPT_DIR}/.."
BUILD_DIR="${ROOT_DIR}/build/macos-cmake"

cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
cmake --build "${BUILD_DIR}" --config Release

echo "Built kbeyond~ into ${BUILD_DIR}"
