#!/bin/bash

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
DEV_DIR="$SCRIPT_DIR/.."

IS_COMMIT_HOOK=0
STRICT_MODE=1

source "$DEV_DIR"/internal/checks-main.sh
