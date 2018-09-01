#!/bin/bash

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

IS_COMMIT_HOOK=0
STRICT_MODE=1

source "$SCRIPT_DIR"/checks-main.sh
