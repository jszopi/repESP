#!/usr/bin/env python3

from repESP import resp, resp_helpers, rep_esp, charges
from repESP.field_comparison import difference, rms_and_rep

import resp_parser

import argparse
import os
import sys
import shutil


help_description = """
    Investigate the dependence of the ESP fit on the charge on one or two atoms
    """

parser = argparse.ArgumentParser(
    parents=[resp_parser.parser],
    description=help_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("esp_file",
                    help=resp_parser.esp_file_help,
                    metavar="FILENAME")

parser.add_argument("--monitor",
                    help="""labels of atoms which charges are to be monitored
                    while the charge on atom1 (and optionally atom2) are varied""",
                    type=int, nargs="*", metavar="LABELS", default=[])

atom1_group = parser.add_argument_group(
    title="options regarding the first varied atom",
    description="Charge on this atom will be varied"
)

atom1_group.add_argument(
    "atom1",
    help="""label of the first atom which charge is to be varied""",
    type=int,
    metavar="LABEL"
)

atom1_group.add_argument(
    "--equivalent1",
    help="""If in the molecule there are atoms equivalent to the atom1, specify
    their labels. This ensures that the atoms cannot vary independently.""",
    type=int,
    nargs="*",
    metavar="LABELS",
    default=[]
)

atom1_group.add_argument(
    "--limits1",
    help="""range of atom1 charge values to be sampled""",
    nargs=2,
    type=float,
    default=(-1, 1),
    metavar=("LOWER", "UPPER")
)

atom1_group.add_argument(
    "--sampling1",
    help="""number of data points to be sampled for atom1 charges""",
    type=float,
    default=11,
    metavar="POINTS"
)

atom2_group = parser.add_argument_group(
    title="options regarding the second varied atom",
    description="Optionally the charge on another atom can be simultanously varied."
)

atom2_group.add_argument(
    "--atom2",
    help="""label of the second atom which charge is to be varied""",
    type=int,
    metavar="LABEL"
)

atom2_group.add_argument(
    "--equivalent2",
    help="""If in the molecule there are atoms equivalent to the atom1, specify
        their labels. This ensures that the atoms cannot vary independently.""",
    type=int,
    nargs="*",
    metavar="LABELS",
    default=[]
)

atom2_group.add_argument(
    "--limits2",
    help="""range of atom2 charge values to be sampled""",
    nargs=2,
    type=float,
    default=(-1, 1),
    metavar=("LOWER", "UPPER")
)

atom2_group.add_argument(
    "--sampling2",
    help="""number of data points to be sampled for atom2 charges""",
    type=float,
    default=11,
    metavar="POINTS"
)

args = parser.parse_args()
