#!/usr/bin/env python3

from repESP import resp, resp_helpers, rep_esp, charges
from repESP.field_comparison import difference, rms_and_rep
from repESP.resp import get_atom_signature

from numpy import linspace, meshgrid

import resp_parser

import argparse
import csv
import itertools
import os
import shutil
import sys


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

parser.add_argument("-o", "--output",
                    help="file to write the output in the csv format",
                    type=str, metavar="FILENAME")

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

def interpret(molecule, charge_dict, vary_label1, vary_label2=None):
    if vary_label2 is None:
        dictio = charge_dict(1)  # Example number to get the dict
    else:
        dictio = charge_dict(1, 2)  # Example numbers to get the dict
    print("\nCharges on these atoms will be varied:")
    for vary_label in vary_label1, vary_label2:
        if vary_label is None:
            break
        print('*', molecule[vary_label-1])
        equiv = [label for label in dictio if
                 dictio[label] == dictio[vary_label] and label != vary_label]
        if equiv:
            print("  with the following atoms equivalenced to it:")
            for equiv_label in sorted(equiv):
                print("  -", molecule[equiv_label-1])
    print("\nSee below for equivalence information of other atoms.")


def get_monitored(molecule, labels):
    return [atom.charges['resp'] for i, atom in enumerate(molecule) if i + 1 in labels]


def sample_charges(all_charge_variations, inp_charges_func, temp_dir_func):

    result = []
    print()

    for i, varied_charges in enumerate(all_charge_variations):

        inp_charges = inp_charges_func(varied_charges)
        temp_dir = temp_dir_func(varied_charges)

        _molecule = resp.run_resp(
            args.respin_location,
            temp_dir,
            resp_type='dict',
            inp_charges=inp_charges,
            esp_fn=args.esp_file,
            check_ivary=i==0  # Only for the first iteration
        )

        rms, rrms, _ = rms_and_rep(info_from_esp.field, _molecule, 'resp')

        sys.stdout.write("\rSampling progress: {0:.2f} %".format(
            100*(i+1)/len(all_charge_variations)))
        sys.stdout.flush()

        result.append([
            rms,
            rrms,
            *get_monitored(_molecule, args.monitor)
        ])

        # For large resolutions this directory can become impractically large
        shutil.rmtree(temp_dir)

    print()
    return result


def get_csv_header(monitored, atom1, atom2=None):
    get_charge_header = lambda label: "Charge on {}".format(label)
    header_common = ["RMS", "RRMS", *list(map(get_charge_header, monitored))]
    if atom2 is None:
        return [get_charge_header(atom1)] + header_common
    else:
        return [get_charge_header(atom1), get_charge_header(atom2)] + header_common


def one_charge_variation():

    charge_dict = lambda x: {a: x for a in args.equivalent1 + [args.atom1]}

    interpret(info_from_esp.molecule, charge_dict, args.atom1)

    all_charge_variations = linspace(args.limits1[0], args.limits1[1], num=args.sampling1)

    inp_charges_func = lambda varied_charge: resp.charges_from_dict(
        charge_dict(varied_charge),
        len(info_from_esp.molecule)
    )

    temp_dir_func = lambda varied_charge: temp_dir + "{0}{1:+.3f}".format(
        get_atom_signature(info_from_esp.molecule, args.atom1), varied_charge,
    )

    results = sample_charges(all_charge_variations, inp_charges_func, temp_dir_func)

    charges_with_results = [[charges, *result] for charges, result in zip(all_charge_variations, results)]
    csv_header = get_csv_header(args.monitor, args.atom1)

    return csv_header, charges_with_results


def two_charge_variation():

    charge_dict = lambda x, y: {
        **{a: x for a in args.equivalent1 + [args.atom1]},
        **{b: y for b in args.equivalent2 + [args.atom2]}
    }

    interpret(info_from_esp.molecule, charge_dict, args.atom1, args.atom2)

    all_charge_variations = list(itertools.product(
        linspace(args.limits1[0], args.limits1[1], num=args.sampling1),
        linspace(args.limits2[0], args.limits2[1], num=args.sampling2)
    ))

    inp_charges_func = lambda varied_charges: resp.charges_from_dict(
        charge_dict(*varied_charges),
        len(info_from_esp.molecule)
    )

    temp_dir_func = lambda varied_charges: temp_dir + "{0}{1:+.3f}-{2}{3:+.3f}".format(
        get_atom_signature(info_from_esp.molecule, args.atom1), varied_charges[0],
        get_atom_signature(info_from_esp.molecule, args.atom2), varied_charges[1],
    )

    results = sample_charges(all_charge_variations, inp_charges_func, temp_dir_func)

    charges_with_results = [[*charges, *result] for charges, result in zip(all_charge_variations, results)]
    csv_header = get_csv_header(args.monitor, args.atom1, args.atom2)

    return csv_header, charges_with_results


args = parser.parse_args()

input_esp = args.respin_location + "/" + args.esp_file

temp_dir = "fit-dependence_temp_dir/"

if args.output and os.path.exists(args.output):
    raise FileExistsError("Output file exists: " + args.output)

if os.path.exists(temp_dir):
    shutil.rmtree(temp_dir)

os.mkdir(temp_dir)

# Read the .esp file
info_from_esp = resp_helpers.G09_esp(input_esp)
# Write the .esp file in the correct format expected by the `resp` program
info_from_esp.field.write_to_file(temp_dir + "corrected.esp", info_from_esp.molecule)

csv_header, charges_with_results = one_charge_variation() if args.atom2 is None else two_charge_variation()

if args.output is None:
    print("\nESP fit dependence scan results:")
    print(",".join("{:>13}".format(x) for x in csv_header))
    for result_line in charges_with_results:
        print(",".join("{:>13.8f}".format(x) for x in result_line))
else:
    with open(args.output, "w") as out:
        csv_writer = csv.writer(out)
        csv_writer.writerow(csv_header)
        csv_writer.writerows(charges_with_results)

shutil.rmtree(temp_dir)
