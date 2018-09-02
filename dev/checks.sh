#!/bin/bash

DEV_DIR="$(dirname "$(readlink -f "$0")")"

IS_COMMIT_HOOK=0
STRICT_MODE=0

source "$DEV_DIR"/internal/checks-main.sh
