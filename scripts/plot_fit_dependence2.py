#!/usr/bin/env python3

from plot_fit_dependence_common import get_parser

if __name__ == "__main__":

    parser, plot_appearance_group = get_parser(isTwoAtoms=True)

    args = parser.parse_args()
