import ipdb
from cube_helpers import Atom, Cube

esp_type_in_log = {
    ' Merz-Kollman atomic radii used.': 'mk',
    ' Francl (CHELP) atomic radii used.': 'chelp',
    ' Breneman (CHELPG) radii used.': 'chelpg',
    }

esp_charges = esp_type_in_log.values()


# First 8 characters of the line following charge output in various input files
charge_termination_line = {
    'sumviz': '--------',
    'log': ' Sum of ',
    }


class InputFortmatError(Exception):
    pass


def get_charges(charge_type, filename, molecule):
    """Update the molecule with charges

    Only charges calculated directly by Gaussian are currently supported. The
    charges should be given in a Gaussian output file (.log or .out). In the
    future, checkpoint and formatted checkpoint formats will be supported.
    """
    if filename[-4:] in ['.log', '.out']:
        _charges_from_log(charge_type, filename, molecule)
    elif filename[-4:] in ['.chk', '.fchk']:
        raise NotImplementedError('File extension {0} currently not supported.'
                                  .format(filename[-4]))
    elif filename[-7:] == '.sumviz' and charge_type == 'aim':
        _charges_from_sumviz(filename, molecule)
    else:
        raise NotImplementedError('File extension {0} is not supported.'
                                  .format(filename[-4]))


def _charges_from_log(charge_type, filename, molecule):
    """Update the molecule with charges from Gaussian output

    Note that if the output file contains information about charges in more
    than one place, only the last one will be used. Also, the atom list is
    assumed to be in order.
    """
    with open(filename, 'r') as f:
        _goto_occurence_in_log(charge_type, f, -1)
        charges = _get_charges_from_lines(f, 'log', molecule)

        for atom, charge in zip(molecule, charges):
            atom.charges[charge_type] = float(charge)


def _goto_occurence_in_log(charge_type, file_object, occurence):
    """Go to the selected occurence of input about charges in a log file.

    Occurence is the index to a list containing all occurences of the given
    charge type, so should be 0 for the first occurence and -1 for the last.
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
        if line == _charge_section_header_in_log(charge_type):
            result.append(offset)
        # The information about the type of ESP charges is gathered separately
        if charge_type in esp_charges and line in esp_type_in_log:
            esp_types.append(esp_type_in_log[line])

    if charge_type in esp_charges:
        # Verify if all ESP charge output has been recognized correctly
        if len(esp_types) != len(result):
            raise InputFortmatError('Information about the type of some '
                                    'ESP charges was not recognized.')
        # Filter only the requested ESP charge type
        result = [elem for i, elem in enumerate(result) if
                  esp_types[i] == charge_type]

    if not result:
        raise InputFortmatError("Output about charge type '{0}' not found."
                                .format(charge_type))

    try:
        file_object.seek(result[occurence])
    except IndexError:
        raise IndexError(
            "Cannot find occurence '{0}' in a list of recognized pieces of "
            "output about charges, whose length is {1}.".format(occurence,
                                                                len(result)))

    # Skip an unnecessary line
    file_object.readline()


def _charge_section_header_in_log(charge_type):
    if charge_type == 'mulliken':
        name = 'Mulliken'
    elif charge_type in esp_charges:
        name = 'ESP'
    else:
        raise NotImplementedError("Charge type '{0}' is not implemented."
                                  .format(charge_type))
    return ' ' + name + ' charges:'


def _charges_from_sumviz(filename, molecule):
    """Update the molecule with charges from AIMAll output."""
    with open(filename, 'r') as f:
        while f.readline().rstrip('\n') != 'Some Atomic Properties:':
            pass
        # Skip irrelevant lines
        for i in range(9):
            f.readline()

        charges = _get_charges_from_lines(f, 'sumviz', molecule)
        for atom, charge in zip(molecule, charges):
            atom.charges['aim'] = float(charge)


def _log_charge_line(line):
    return line.split()


def _sumviz_charge_line(line):
    letter_and_label, charge, *other = line.split()
    # These should actually be a regex for letters and numbers
    letter = letter_and_label[0]
    label = letter_and_label[1]
    return label, letter, charge


def _get_charges_from_lines(file_object, input_type, molecule):
    charges = []
    for i, atom in enumerate(molecule):
        label, letter, charge = globals()['_' + input_type + '_charge_line'](
            file_object.readline())
        # Check if the atom identities agree between atom list and input
        if letter != atom.identity:
            raise InputFortmatError(
                'Atom {0} in atom list is given as {1} but input file '
                'expected {2}'.format(int(label)+1, atom.identity, letter))
        charges.append(charge)

    # Check if the atom list terminates after as many atoms as expected
    next_line = file_object.readline()
    if next_line[:8] != charge_termination_line[input_type]:
        raise InputFortmatError('Expected end of charges ( \'----------\')'
                                ', instead got: ' + next_line)

    return charges
