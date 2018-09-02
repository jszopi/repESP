#!/bin/bash

DEV_DIR="$(dirname "$(readlink -f "$0")")"
source "$DEV_DIR"/internal/common.sh

run_coverage
TESTS_RC="$?"

if [[ "$TESTS_RC" -ne 0 ]] ; then
    exit "$TESTS_RC"
fi

echo -e "\nTest coverage compared to snapshot:"

COVERAGE_SNAPSHOT="$(run_coverage_snapshot)"
NON_IDENTICAL_LINES="$(echo "$COVERAGE_SNAPSHOT" | grep -v "^\[ == ]")"

if [[ -z "$NON_IDENTICAL_LINES" ]] ; then
    echo "Identical as snapshot."
else
    echo "$COVERAGE_SNAPSHOT"
fi

echo -e "\nSee \`cover/index.html\` for a html version including line breakdown."
