SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
MYPY_SNAPSHOT="$SCRIPT_DIR/mypy_snapshot.out"

replace_numbers() {
    sed -r 's/[0-9]+/0/g'
}

run_mypy() {
    pushd "$SCRIPT_DIR"/.. > /dev/null
    mypy repESP/*.py
    # TODO: mypy check on the tests doesn't seem to be working:
    # mypy tests/*.py
    popd > /dev/null
}

run_tests() {
    nosetests3 -s "$SCRIPT_DIR"/.. 2>&1
}
