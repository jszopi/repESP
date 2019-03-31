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
    mypy repESP
    mypy tests
    popd > /dev/null
}

run_tests() {
    nosetests3 --rednose -s "$ROOT_DIR" 2>&1
}

run_coverage() {
    nosetests3 -s --with-coverage --cover-package repESP --cover-html 2>&1
}

get_coverage() {
    run_coverage  | grep "^repESP" | tr -s " " | cut -f1,4 -d" "
}

run_coverage_snapshot() {
    "$INTERNAL_DIR"/coverage_snapshot.py "$COVERAGE_SNAPSHOT" <(get_coverage)
}
