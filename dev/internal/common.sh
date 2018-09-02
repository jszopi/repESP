if [[ -z "$DEV_DIR" ]] ; then
    echo "Calling script did not specify DEV_DIR. Aborting."
    exit 1
fi

ROOT_DIR="$DEV_DIR/.."
INTERNAL_DIR="$DEV_DIR/internal"

MYPY_SNAPSHOT="$INTERNAL_DIR/mypy_snapshot.out"
COVERAGE_SNAPSHOT="$INTERNAL_DIR/coverage_snapshot.out"

replace_numbers() {
    sed -r 's/[0-9]+/0/g'
}

run_mypy() {
    pushd "$ROOT_DIR" > /dev/null
    mypy repESP/*.py
    # TODO: mypy check on the tests doesn't seem to be working:
    # mypy tests/*.py
    popd > /dev/null
}

run_tests() {
    nosetests3 -s "$ROOT_DIR" 2>&1
}
