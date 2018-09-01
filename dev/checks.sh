#!/bin/bash

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
source "$SCRIPT_DIR"/common.sh

if [[ "$1" == "dev" ]] ; then
    DEV_MODE=1
    echo "Running checks in dev mode..."
else
    DEV_MODE=0
    echo "Running checks..."
fi

echo

# UNIT TESTS

TESTS_OUTPUT="$(run_tests)"
TESTS_RC="$?"

if [[ "$TESTS_RC" -ne 0 ]] ; then
    echo "Unit tests failed with return code $TESTS_RC and output:"
    echo "$TESTS_OUTPUT"
    if [[ "$DEV_MODE" -eq 0 ]] ; then
        exit 1
    fi
else
    echo "Unit tests OK"
fi

echo

# STATIC ANALYSIS

## MYPY SNAPSHOT

if [[ "$DEV_MODE" -eq 0 ]] ; then
    MYPY_DIFF=$(diff -uN "$MYPY_SNAPSHOT" <(run_mypy))
else
    MYPY_DIFF=$(diff -uN <(cat "$MYPY_SNAPSHOT" | replace_numbers) <(run_mypy | replace_numbers))
fi

if [[ -n "$MYPY_DIFF" ]] ; then
    echo -e "mypy snapshot check failed with the following diff:\n"
    echo "$MYPY_DIFF"
    if [[ "$DEV_MODE" -eq 0 ]] ; then
        exit 2
    else
        echo -e "\n(run without dev mode to see line numbers)"
    fi
else
    echo "mypy snapshot check OK"
fi

echo
