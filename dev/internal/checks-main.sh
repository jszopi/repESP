if [[ -z "$DEV_DIR" ]] ; then
    echo "Calling script did not specify DEV_DIR. Aborting."
    exit 1
fi

source "$DEV_DIR"/internal/common.sh

ALL_TESTS_PASSED=1

if [[ "$STRICT_MODE" -eq 1 ]] ; then
    echo "Running checks in strict mode..."
else
    echo "Running checks in dev mode..."
fi

# UNIT TESTS

TESTS_OUTPUT="$(run_tests)"
TESTS_RC="$?"

if [[ "$TESTS_RC" -ne 0 ]] ; then
    echo "[FAIL] Unit tests: failed with return code $TESTS_RC and output:"
    echo "$TESTS_OUTPUT"
    if [[ "$STRICT_MODE" -eq 1 ]] ; then
        exit 2
    fi

    echo
    ALL_TESTS_PASSED=0
else
    echo "[PASS] Unit tests"

    ## TEST COVERAGE SNAPSHOT

    COVERAGE_DIFF=$(run_coverage_snapshot)
    COVERAGE_RC="$?"

    if [[ "$COVERAGE_RC" -ne 0 ]] ; then
        echo "[FAIL] Test coverage: failed with return code $COVERAGE_RC and output:"
        echo "$COVERAGE_DIFF"
        if [[ "$STRICT_MODE" -eq 1 ]] ; then
            exit 3
        fi

        echo
        ALL_TESTS_PASSED=0
    else
        echo "[PASS] Test coverage"
    fi

fi

# STATIC ANALYSIS

## MYPY SNAPSHOT

if [[ "$STRICT_MODE" -eq 1 ]] ; then
    MYPY_DIFF=$(diff -uN "$MYPY_SNAPSHOT" <(run_mypy))
else
    MYPY_DIFF=$(diff -uN <(cat "$MYPY_SNAPSHOT" | replace_numbers) <(run_mypy | replace_numbers))
fi

if [[ -n "$MYPY_DIFF" ]] ; then
    echo -e "[FAIL] mypy snapshot check: failed with the following diff:\n"
    echo "$MYPY_DIFF"

    if [[ "$STRICT_MODE" -eq 1 ]] ; then
        exit 4
    else
        echo -e "\n(run without dev mode to see line numbers)\n"
    fi
    ALL_TESTS_PASSED=0
else
    echo "[PASS] mypy snapshot check"
fi

if [[ "$IS_COMMIT_HOOK" -eq 1 ]] && [[ "$ALL_TESTS_PASSED" -eq 0 ]] ; then
    echo
    exec < /dev/tty
    read -p "Press any key to proceed with the commit, Ctrl+c to abort."
    exec <&-
    echo
fi
