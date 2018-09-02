#!/bin/bash

DEV_DIR="$(dirname "$(readlink -f "$0")")"
source "$DEV_DIR"/internal/common.sh

run_mypy > "$MYPY_SNAPSHOT"
