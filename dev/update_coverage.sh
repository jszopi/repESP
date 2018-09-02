#!/bin/bash

DEV_DIR="$(dirname "$(readlink -f "$0")")"
source "$DEV_DIR"/internal/common.sh

COVERAGE_OUTPUT=$(run_coverage_snapshot)

get_coverage > "$COVERAGE_SNAPSHOT"
