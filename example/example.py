#!/usr/bin/env python3

# This script is an example use case for the library. See html documentation
# for a description (or repESP/__init__.py file if you're unable to generate
# the documentation at the moment).

# 1. Extract ESP field at CHelpG fitting points:

from repESP.esp_util import parse_gaussian_esp

with open("../data/NMe3H_plus/NMe3H_plus_chelpg.esp") as f:
    gaussian_esp_data = parse_gaussian_esp(f)

molecule = gaussian_esp_data.molecule
esp_field = gaussian_esp_data.field

print("Extracted ESP field at CHelpG fitting points (not shown).")
print("=========================================================")
print()

print("Extracted molecule from ESP file (charges are CHelpG):")
print("======================================================")
for atom in molecule.atoms:
    print(atom)
print()

# 2. Parse NPA charges from Gaussian output:

from repESP.gaussian_format import get_charges_from_log, NpaChargeSectionParser

with open("../data/NMe3H_plus/NMe3H_plus_nbo.log") as f:
    npa_charges = get_charges_from_log(f, NpaChargeSectionParser(), molecule)

print("Extracted NPA charges:")
print("======================")
for charge in npa_charges:
    print(charge)
print()

# 3. Extract equivalence information from respin files:

from repESP.respin_format import parse_respin, get_equivalence_from_two_stage_resp_ivary

with open("../data/NMe3H_plus/NMe3H_plus_mk.respin1") as respin1_file:
    respin1 = parse_respin(respin1_file)

with open("../data/NMe3H_plus/NMe3H_plus_mk.respin2") as respin2_file:
    respin2 = parse_respin(respin2_file)

equivalence = get_equivalence_from_two_stage_resp_ivary(respin1.ivary, respin2.ivary)

print("Extracted equivalence information:")
print("==================================")
print(equivalence.describe(molecule))

# 4. Fit hydrogens using NPA as initial charges:

from repESP.esp_util import EspData
from repESP.resp_wrapper import fit_hydrogens_only

esp_data = EspData.from_gaussian(gaussian_esp_data)

my_charges = fit_hydrogens_only(
    esp_data,
    equivalence,
    molecule,
    total_charge=1,  # cation, +1 total charge
    initial_charges=npa_charges  # will be used for non-hydrogen atoms
)

print("My new charges were computed:")
print("=============================")
for charge in my_charges:
    print(charge)
print()

# 5. Reproduce ESP and calculate fit statistics:

from repESP.charges import AtomWithCoordsAndCharge
from repESP.types import Molecule

molecule_with_my_charges = Molecule([
    AtomWithCoordsAndCharge(
        atom.atomic_number,
        atom.coords,
        my_charge
    ) for atom, my_charge in zip(molecule.atoms, my_charges)
])

from repESP.calc_fields import esp_from_charges, calc_rms_error, calc_relative_rms_error

reproduced_esp = esp_from_charges(esp_field.mesh, molecule_with_my_charges)

print("Reproduced ESP field from my new charges (not shown).")
print("=====================================================")
print()

rms = calc_rms_error(esp_field.values, reproduced_esp.values)
rrms = calc_relative_rms_error(esp_field.values, reproduced_esp.values)

print("Calculated fit statistics w.r.t. actual ESP field.")
print("==================================================")
print(" RMS: {:.6f} a.u.".format(rms))
print("RRMS: {:.6f}".format(rrms))
