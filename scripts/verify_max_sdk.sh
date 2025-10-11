#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: ./scripts/verify_max_sdk.sh [path]

Checks that the provided Max SDK path contains the headers required by kbeyond~.
If no path is supplied, MAX_SDK_ROOT is used. The script prints diagnostics to
help troubleshoot common configuration mistakes (e.g. forgotten unzip of the SDK).
USAGE
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
    usage
    exit 0
fi

SDK_PATH="${1:-${MAX_SDK_ROOT:-}}"

if [[ -z "${SDK_PATH}" ]]; then
    echo "ERROR: no Max SDK path supplied. Pass it as an argument or set MAX_SDK_ROOT." >&2
    exit 1
fi

SDK_PATH="$(cd "${SDK_PATH}" 2>/dev/null && pwd)" || {
    echo "ERROR: \"${SDK_PATH}\" does not exist or cannot be accessed." >&2
    exit 1
}

MAX_INCLUDES="${SDK_PATH}/source/c74support/max-includes"
MSP_INCLUDES="${SDK_PATH}/source/c74support/msp-includes"

missing=0
if [[ ! -f "${MAX_INCLUDES}/ext.h" ]]; then
    echo "ERROR: Missing ext.h at ${MAX_INCLUDES}." >&2
    missing=1
fi

if [[ ! -f "${MSP_INCLUDES}/z_dsp.h" ]]; then
    echo "ERROR: Missing z_dsp.h at ${MSP_INCLUDES}." >&2
    missing=1
fi

if [[ ${missing} -ne 0 ]]; then
    cat >&2 <<'HINT'
The path exists but appears incomplete. Make sure you downloaded and unzipped the
full Max SDK (available from https://cycling74.com/downloads/sdk ). The correct
layout includes source/c74support/max-includes and msp-includes.
HINT
    exit 1
fi

echo "Max SDK detected at: ${SDK_PATH}"
echo "  max-includes : ${MAX_INCLUDES}"
echo "  msp-includes : ${MSP_INCLUDES}"
echo "All required headers are present."
