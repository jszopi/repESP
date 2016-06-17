from .cube_helpers import Atom, Cube, InputFormatError

import re
import math
from fortranformat import FortranRecordWriter

esp_type_in_log = {
    ' Merz-Kollman atomic radii used.': 'mk',
    ' Francl (CHELP) atomic radii used.': 'chelp',
    ' Breneman (CHELPG) radii used.': 'chelpg',
    ' Generate Potential Derived Charges using the Hu-Lu-Yang model:': 'hly',
    }

esp_charges = esp_type_in_log.values()


def get_rms_from_log(filename):
    # Simple regex for floating point numbers will suffice, as they are
    # non-negative and in decimal format:
    # http://stackoverflow.com/a/4703409
    rms_re = re.compile(" Charges from ESP fit, RMS=\s+(\d+\.\d+) RRMS=\s+"
                        "(\d+\.\d+):$")
    found_line = None
    with open(filename, 'r') as file_object:
        for line_num, line in enumerate(file_object):
            if rms_re.match(line) is not None:
                if found_line is not None:
                    raise InputFormatError(
                        "Multiple lines match the ESP summary pattern: lines "
                        "{0} and {1}.".format(found_line+1, line_num+1))
                rms_line = rms_re.match(line)
                found_line = line_num
        if rms_line is None:
            raise InputFormatError("No ESP fit summary found.")
        return float(rms_line.group(1)), float(rms_line.group(2))


def dump_charges_to_qout(molecule, charge_type, filename):
    line = FortranRecordWriter("8F10.6")
    charges = [atom.charges[charge_type] for atom in molecule]
    with open(filename, 'w') as f:
        print(line.write(charges), file=f)


def update_with_charges(charge_type, filename, molecule, verbose=True):
    """Update the molecule with charges

    Only charges calculated directly by Gaussian are currently supported. The
    charges should be given in a Gaussian output file (.log or .out). In the
    future, checkpoint and formatted checkpoint formats will be supported.

    Note that if Gaussian output file contains information about charges in
    more than one place, only the last one will be used. Also, the atom list is
    assumed to be in order.
    """
    if filename[-4:] in ['.log', '.out']:
        _get_charges(charge_type, filename, 'log', molecule)
    elif filename[-7:] == '.sumviz' and charge_type == 'aim':
        _get_charges('aim', filename, 'sumviz', molecule)
    elif filename[-4:] == '.dat' and charge_type == 'aim':
        if verbose:
            print("WARNING: The QTAIM charges obtained from Henkelman group's "
                  "`bader` program differed from those from `AIMAll` for "
                  "methane by more than 0.01e, even with a very fine cube "
                  "grid. Extracting charges from `.dat` files is hence "
                  "currently not recommended.")
        _get_charges('aim', filename, 'dat', molecule)
    elif filename[-5:] == ".qout":
        if verbose:
            print("WARNING: qout file doesn't contain information about the "
                  "identities of atoms. You should check if it has been "
                  "generated (likely by the `resp` program) using the same "
                  "order of atoms as that in the molecule being updated.")
        _get_charges(charge_type, filename, 'qout', molecule)
    elif filename[-4:] in ['.chk', '.fchk']:
        raise NotImplementedError('File extension {0} currently not supported.'
                                  .format(filename[-4]))
    else:
        raise NotImplementedError('File extension {0} is not supported.'
                                  .format(filename[-4]))


def compare_charges(charge_type1, charge_type2, molecule, molecule2=None,
                    thresh=0.05):
    """Check if two types of charges differ by more than a threshhold

    The charges can be in the same molecule or another molecule object with the
    same atoms.
    """
    output = []
    if molecule2 is None:
        molecule2 = molecule
    for atom1, atom2 in zip(molecule, molecule2):
        if atom1 != atom2:
            raise ValueError("The given molecules contain different atoms!\n"
                             "{0}\n{1}".format(atom1, atom1))
        # May rise KeyError if any atom doesn't have the right charges
        diff = atom1.charges[charge_type1] - atom2.charges[charge_type2]
        # The second condition tests for equality if the threshold is not zero
        if abs(diff) > thresh or (thresh and math.isclose(abs(diff), thresh)):
            output.append("{0} differ by {1: .3f}".format(atom1, diff))
    return '\n'.join(output)


def _get_charges(charge_type, filename, input_type, molecule):
    """Update the molecule with charges."""
    with open(filename, 'r') as file_object:
        if input_type == 'qout':
            charges = _get_charges_from_qout(file_object)
        else:
            globals()['_goto_in_' + input_type](charge_type, file_object)
            charges = _get_charges_from_lines(charge_type, file_object,
                                              input_type, molecule)
        _update_molecule_with_charges(molecule, charges, charge_type)


def _get_charges_from_qout(file_object):
    charges = []
    for line in file_object:
        charges += [float(elem) for elem in line.split()]
    return charges


def _goto_in_log(charge_type, file_object, occurrence=-1):
    """Go to the selected occurrence of input about charges in a log file.

    Occurrence is the index to a list containing all occurrences of the given
    charge type, so should be 0 for the first occurrence and -1 for the last.
    Code based on: http://stackoverflow.com/a/620492
    """
    offset = 0
    result = []
    esp_types = []

    for line in file_object:
        offset += len(line)
        line = line.rstrip('\n')
        # All ESP charges are added here, as they cannot be distinguished just
        # by the header
        if _charge_section_header_in_log(charge_type) in line.rstrip():
            result.append(offset)
        # The information about the type of ESP charges is gathered separately
        if charge_type in esp_charges and line in esp_type_in_log:
            esp_types.append(esp_type_in_log[line])

    if charge_type in esp_charges:
        # Verify if all ESP charge output has been recognized correctly
        if len(esp_types) != len(result):
            raise InputFormatError('Information about the type of some '
                                   'ESP charges was not recognized.')
        # Filter only the requested ESP charge type
        result = [elem for i, elem in enumerate(result) if
                  esp_types[i] == charge_type]

    if not result:
        raise InputFormatError("Output about charge type '{0}' not found."
                               .format(charge_type))

    try:
        file_object.seek(result[occurrence])
    except IndexError:
        raise IndexError(
            "Cannot find occurrence '{0}' in a list of recognized pieces of "
            "output about charges, whose length is {1}.".format(occurrence,
                                                                len(result)))

    # Skip unnecessary lines
    lines_count = 1
    if charge_type == 'nbo':
        lines_count = 5
    if charge_type in esp_type_in_log.values():
        lines_count = 2
    for counter in range(lines_count):
        file_object.readline()


def _goto_in_sumviz(charge_type, file_object):
    while file_object.readline().rstrip('\n') != 'Some Atomic Properties:':
        pass
    # Skip irrelevant lines
    for i in range(9):
        file_object.readline()


def _goto_in_dat(charge_type, file_object):
    for i in range(2):
        file_object.readline()


def _charge_section_header_in_log(charge_type):
    if charge_type == 'mulliken':
        return ' Mulliken charges:'
    elif charge_type in esp_charges:
        return ' Charges from ESP fit,'
    elif charge_type == 'nbo':
        return ' Summary of Natural Population Analysis:'
    else:
        raise NotImplementedError("Charge type '{0}' is not implemented."
                                  .format(charge_type))


def _charge_termination_line(input_type, charge_type):
    """Returns the first 8 characters of the charge section termination line"""
    if input_type == 'log':
        if charge_type == 'nbo':
            return [' =======']
        else:
            # The latter happens when IOp(6/50=1) is requested
            return [' Sum of ', ' Charges']
    elif input_type == 'sumviz':
        return ['--------']
    elif input_type == 'dat':
        return [' -------']
    else:
        raise NotImplementedError("Combination of input file type '{0}' and "
                                  "charge type '{1}' is not implemented."
                                  .format(input_type, charge_type))


def _update_molecule_with_charges(molecule, charges, charge_type):
    if len(molecule) != len(charges):
        raise ValueError("The number of charges provided ({0}) is not equal to"
                         " the number of atoms in molecule ({1}).".format(
                             len(charges), len(molecule)))
    for atom, charge in zip(molecule, charges):
        atom.charges[charge_type] = charge


def _log_charge_line(line, charge_type):
    if charge_type == 'nbo':
        letter, label, charge, *other = line.split()
    else:
        label, letter, charge = line.split()
    label = int(label)
    charge = float(charge)
    return label, letter, charge


def _sumviz_charge_line(line, charge_type):
    letter_and_label, charge, *other = line.split()
    # These should actually be a regex for letters and numbers
    letter = letter_and_label[0]
    label = int(letter_and_label[1])
    charge = float(charge)
    return label, letter, charge


def _dat_charge_line(line, charge_type):
    label, x, y, z, charge, *other = line.split()
    label = int(label)
    charge = float(charge)
    return label, None, charge


def _get_charges_from_lines(charge_type, file_object, input_type, molecule):
    """Extract charges from the charges section in output

    Parameters
    ----------
    file_object : File
        The file from which the charges are to be extracted. The file is
        expected to be set to the position of the start of charges section,
        e.g. with the _goto_in_log helper.
    input_type : str
        Currently implemented is reading lines from Gaussian ('log') and AIM
        ('sumviz') output files.
    molecule : Molecule
        The molecule to which the charges relate. Note that the molecule will
        not be updated with the charges, this must be done separately by the
        caller.

    Returns
    -------
    List[float]
        List of charges in order of occurrence in output file.

    Raises
    ------
    NotImplementedError
        Raised when an unsupported input file type is requested.
    InputFormatError
        Raised when the order of atoms is not as expected from the Molecule or
        the length of the charges section is different than expected.

    """
    charges = []
    for i, atom in enumerate(molecule):
        try:
            # Input type-specific extraction performed by specialist function
            label, letter, charge = globals()[
                '_' + input_type + '_charge_line'](file_object.readline(),
                                                   charge_type)
        except KeyError:
            raise NotImplementedError(
                "Reading charges from an input file of type '{0} 'is not "
                "supported.".format(input_type))

        # Check if the labels are in order
        if label is not None and label != i + 1:
            raise InputFormatError(
                "Charge section is not given in order of Gaussian labels. This"
                " may be a feature of the program which generated the charges "
                "output but is not supported in this program.")
        # Check if the atom identities agree between atom list and input
        if letter is not None and letter != atom.identity:
            raise InputFormatError(
                'Atom {0} in atom list is given as {1} but input file '
                'expected {2}'.format(int(label)+1, atom.identity, letter))

        if charge_type == 'aim' and file_object.name[-4:] == '.dat':
            charge = atom.atomic_no - charge

        charges.append(charge)

    # Check if the atom list terminates after as many atoms as expected from
    # the Molecule object given
    next_line = file_object.readline()
    expected = _charge_termination_line(input_type, charge_type)
    if next_line[:8] not in expected:
        expected = "' or '".join(expected)
        raise InputFormatError(
            "Expected end of charges ('{0}'), instead got: '{1}'".format(
                expected, next_line[:8]))

    return charges
