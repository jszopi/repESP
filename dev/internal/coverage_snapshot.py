#!/usr/bin/env python3

import argparse
import math
from typing import List, Dict, NamedTuple, Optional, TextIO
import sys


CoveredFiles = Dict[str, float]


class FileDiff(NamedTuple):
    filename: str
    old_coverage: Optional[float]
    new_coverage: Optional[float]


def parse_report(file: TextIO) -> CoveredFiles:
    covered_files: CoveredFiles = {}
    for i, line in enumerate(file):
        try:
            line = line.strip()
            assert line[-1] == "%"

            line_split = line[:-1].split()
            assert len(line_split) == 2

            filename, coverage = line_split

            assert filename not in covered_files
            covered_files[filename] = float(coverage)
        except (AssertionError, ValueError):
            print(
                f"Skiping line {i+1} of file {file.name}. The malformed line is:\n{line}",
                file=sys.stderr
            )

    return covered_files


def diff_reports(old_report: CoveredFiles, new_report: CoveredFiles) -> List[FileDiff]:

    all_keys = {*list(old_report.keys()), *list(new_report.keys())}

    return [
        FileDiff(
            filename,
            old_report[filename] if filename in old_report else None,
            new_report[filename] if filename in new_report else None
        ) for filename in sorted(all_keys)
    ]


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


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
            Diff two coverage reports, where each line is "$FILENAME $COVERAGE%".
            The exit code is 1 if the coverage worsened, 0 otherwise.
        """
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
