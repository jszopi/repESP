#!/usr/bin/env python3

from repESP import resp, resp_helpers, rep_esp, charges
from repESP.field_comparison import difference, rms_and_rep

import resp_parser

import argparse
import os
import sys
import shutil


help_description = """
    Calculate the flexibilty limits of an atom evaluated on the given mesh.

    Don't forget to specify the --equivalent option if there are other atoms
    equivalent to the investigated one.

    This script assumes that the ESP fit error as a function of the charge on
    the selected atom has a single minimum and the error increases
    monotonically in both directions away from the minimum. As flexibility is a
    new concept, this may not be the case in all molecules. To study the
    dependence of the ESP fit on the charge of an atom, please use the
    `fit_dependence` script.
    """

parser = argparse.ArgumentParser(
    parents=[resp_parser.parser],
    description=help_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("esp_file",
                    help=resp_parser.esp_file_help,
                    metavar="FILENAME")

parser.add_argument("atom",
                    help="""Label of the atom which flexibility limits are to
                    be evaluated.""",
                    type=int, metavar="LABEL")

parser.add_argument("--equivalent",
                    help="""If in the molecule there are atoms equivalent to
                    the investigated one, specify their labels. This ensures
                    that the atoms cannot vary independently.""",
                    type=int, nargs="*", metavar="LABELS", default=[])

parser.add_argument("--levels",
                    help="""Flexibility levels to be calculated""",
                    type=int, nargs="*", metavar="LEVELS",
                    default=[10])

parser.add_argument("--guess",
                    help="""Guess for the value of the flexibility corresponding
                    to the highest flexibility level. If the specified guess is
                    too small, the optimization may fail, so this number should
                    be the upper limit for the expected flexibility.""",
                    type=float, metavar="VALUE", default=1)

args = parser.parse_args()
input_esp = args.respin_location + "/" + args.esp_file

temp_dir = "flexibility_temp_dir/"

if os.path.exists(temp_dir):
    shutil.rmtree(temp_dir)

os.mkdir(temp_dir)

# Read the .esp file
info_from_esp = resp_helpers.G09_esp(input_esp)
# Write the .esp file in the correct format expected by the `resp` program
info_from_esp.field.write_to_file(temp_dir + "corrected.esp", info_from_esp.molecule)

print("\nRunning unrestrained RESP to fit ESP with equivalence:")
esp_equiv_molecule = resp.run_resp(
    args.respin_location,
    temp_dir + 'unrest',
    resp_type='unrest',
    esp_fn=args.esp_file
)

resp_rms, resp_rrms = rms_and_rep(info_from_esp.field, esp_equiv_molecule, 'resp')[:2]

print("\nThe molecule with equivalenced charges (unrestrained RESP):")

for atom in esp_equiv_molecule:
    atom.print_with_charge('resp')

print()
print(" RMS: {0:.5f}".format(resp_rms))
print("RRMS: {0:.5f}".format(resp_rrms))
print("RMSV: {0:.5f}".format(resp_rms/resp_rrms))

resp_charge = esp_equiv_molecule[args.atom - 1].charges['resp']
upper_charge = resp_charge + args.guess/2
lower_charge = resp_charge - args.guess/2

charge_dict = lambda x: {a: x for a in args.equivalent + [args.atom]}

resp_args = [
    info_from_esp.field,
    args.respin_location,
    temp_dir,
    args.esp_file,
    info_from_esp.molecule,
    args.atom,
    charge_dict
]

upper_rrms, _ = resp.eval_one_charge_resp(upper_charge, *resp_args, False)
lower_rrms, _ = resp.eval_one_charge_resp(lower_charge, *resp_args, False)

flex_limits = []
print("\nFlexibility limits on", info_from_esp.molecule[args.atom - 1])

for level in sorted(args.levels):
    if (1+level/100)*resp_rrms > min(lower_rrms, upper_rrms):
        # Otherwise find_flex would throw an error
        print("Further limits were beyond the sampled range.")
        break

    sol1, sol2, charges1, charges2 = resp.find_flex(
        (1+level/100)*resp_rrms,
        [lower_charge, resp_charge, upper_charge],
        [lower_rrms, resp_rrms, upper_rrms],
        resp_args
    )

    flex_limits.append([sol1, sol2])

    print("{0:>3}% limits: {1: .5f}, {2: .5f}, diff: {3:.5f} ({4:.1f}%)".format(
          level, sol1, sol2, sol2-sol1, 100*abs((sol2-sol1)/resp_charge)
    ))

shutil.rmtree(temp_dir)
