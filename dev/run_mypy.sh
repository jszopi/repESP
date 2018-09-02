#!/bin/bash

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
source "$SCRIPT_DIR"/common.sh

run_mypy
