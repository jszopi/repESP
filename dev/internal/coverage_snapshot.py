#!/usr/bin/env python3.6

import argparse
import math
from typing import Any, List, Dict, NamedTuple, Optional, TextIO
import sys


CoveredFiles = Dict[str, float]


class FileDiff(NamedTuple):
    filename: str
    old_coverage: Optional[float]
    new_coverage: Optional[float]


def parse_report(file: TextIO) -> CoveredFiles:
    covered_files: CoveredFiles = {}
    for line in file:
        line = line.strip()
        assert line[-1] == "%"

        line_split = line[:-1].split()
        assert len(line_split) == 2

        filename, coverage = line_split

        assert filename not in covered_files
        covered_files[filename] = float(coverage)

    return covered_files


def diff_reports(old_report: CoveredFiles, new_report: CoveredFiles) -> List[FileDiff]:

    file_diffs: List[FileDiff] = []
    all_keys = {*list(old_report.keys()), *list(new_report.keys())}

    for filename in sorted(all_keys):
        if filename not in new_report:
            file_diffs.append(FileDiff(filename, old_report[filename], None))
            continue
        if filename not in old_report:
            file_diffs.append(FileDiff(filename, None, new_report[filename]))
            continue
        file_diffs.append(FileDiff(filename, old_report[filename], new_report[filename]))

    return file_diffs


def compare_reports(file_diffs: List[FileDiff]) -> bool:
    new_report_is_worse = False

    for file_diff in file_diffs:
        filename, old_coverage, new_coverage = file_diff

        if new_coverage is None:
            print(f"[DEL ] {filename}: {old_coverage} %")
            new_report_is_worse = True
        elif old_coverage is None:
            print(f"[NEW ] {filename}: {new_coverage} %")
        elif math.isclose(new_coverage, old_coverage):
            print(f"[ == ] {filename}: {old_coverage} %")
        else:
            is_worse = new_coverage < old_coverage
            new_report_is_worse |= is_worse
            tag = "DOWN" if is_worse else " UP "
            print(f"[{tag}] {filename}: {old_coverage} -> {new_coverage} % ")

    return new_report_is_worse


def parse_args() -> List[Any]:
    parser = argparse.ArgumentParser(
        description='Diff two coverage reports, where each line is "$FILENAME $COVERAGE%"'
    )

    parser.add_argument('old_report', metavar='OLD', help="old report", type=argparse.FileType('r'))
    parser.add_argument('new_report', metavar='NEW', help="new report", type=argparse.FileType('r'))

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    old_report = parse_report(args.old_report)
    new_report = parse_report(args.new_report)

    file_diffs = diff_reports(old_report, new_report)
    new_report_is_worse = compare_reports(file_diffs)

    sys.exit(1 if new_report_is_worse else 0)
