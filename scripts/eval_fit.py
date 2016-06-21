#!/usr/bin/env python3

from repESP.field_comparison import rms_and_rep
from repESP import charges
from repESP import resp_helpers

import argparse
import charges_parser
import resp_parser

import os
import shutil

help_description = """Evaluate the fit quality of the given charges on the
    provided mesh of fitting points."""

parser = argparse.ArgumentParser(
    parents=[charges_parser.parser],
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=help_description)

parser.add_argument("esp_filename",
                    help=resp_parser.esp_file_help,
                    metavar="ESP_FILE_NAME")

args = parser.parse_args()

esp_file = resp_helpers.G09_esp(args.esp_filename)
molecule = esp_file.molecule

input_type = charges_parser.input_type(args.input_charge_type)
charges._get_charges(args.input_charge_type, args.input_charge_file,
                     input_type, molecule)

rms, rrms = rms_and_rep(esp_file.field, molecule, args.input_charge_type)[:2]
print(" RMS: {0:.5f}".format(rms))
print("RRMS: {0:.5f}".format(rrms))
print("RMSV: {0:.5f}".format(rms/rrms))
