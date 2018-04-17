#!/usr/bin/env python3

from plot_fit_dependence_common import get_parser, preprocess_args

import pandas

if __name__ == "__main__":

    parser, plot_appearance_group = get_parser(isTwoAtoms=True)

    args = parser.parse_args()
    preprocess_args(args)

    df = pandas.read_csv(args.scan_output)
