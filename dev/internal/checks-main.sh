source "$SCRIPT_DIR"/common.sh

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
        exit 1
    else
        echo
    fi
    ALL_TESTS_PASSED=0
else
    echo "[PASS] Unit tests"
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
        exit 2
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
