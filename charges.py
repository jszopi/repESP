import ipdb
from cube_helpers import Atom, Cube

esp_charges = ['mk', 'chelp', 'chelpg']


class NotImplementedError(Exception):
    pass


class InputFortmatError(Exception):
    pass


def get_charges(charge_type, filename, atoms):
    """Update the atom lists with charges

    Only charges calculated directly by Gaussian are currently supported. The
    charges should be given in a Gaussian output file (.log or .out). In the
    future, checkpoint and formatted checkpoint formats will be supported.
    """
    if filename[-4:] in ['.log', '.out']:
        _charges_from_log(charge_type, filename, atoms)
    elif filename[-4:] in ['.chk', '.fchk']:
        raise NotImplementedError('File extension {0} currently not supported.'
                                  .format(filename[-4]))
    else:
        raise NotImplementedError('File extension {0} is not supported.'
                                  .format(filename[-4]))


def _charges_from_log(charge_type, filename, atoms):
    """Update the atom lists with charges from Gaussian output

    Note that if the output file contains information about charges in more
    than one place, only the last one will be used. Also, the atom list is
    assumed to be in order.
    """
    with open(filename, 'r') as f:
        _goto_occurence_in_log(charge_type, f, -1)
        charges = []

        for i, atom in enumerate(atoms):
            label, letter, charge = f.readline().split()
            # Check if the atom identities agree between atom list and input
            if letter != atom.identity:
                raise InputFortmatError(
                    'Atom {0} in atom list is given as {1} but input file '
                    'expected {2}'.format(int(label)+1, atom.identity, letter))
            charges.append(charge)

        # Check if the atom list terminates after as many atoms as expected
        next_line = f.readline()
        if next_line[:8] != ' Sum of ':
            raise InputFortmatError('Expected end of charges ( \'Sum of ...\')'
                                    ', instead got: ' + next_line)

        for atom, charge in zip(atoms, charges):
            atom.charges[charge_type] = float(charge)


def _goto_occurence_in_log(charge_type, file_object, occurence):
    """Go to the selected occurence of input about charges in a log file.

    Occurence is the index to a list containing all occurences, so should be 0
    for the first occurence and -1 for the last. Code based on:
    http://stackoverflow.com/a/620492
    """
    offset = 0
    result = []

    for line in file_object:
        offset += len(line)
        line = line.rstrip('\n')
        if line == _charge_section_header_in_log(charge_type):
            result.append(offset)

    if result:
        file_object.seek(result[occurence])
        # Skip an unnecessary line
        file_object.readline()
    else:
        raise InputFortmatError("Output about charge type '{0}' not found."
                                .format(charge_type))


def _charge_section_header_in_log(charge_type):
    if charge_type == 'mulliken':
        name = 'Mulliken'
    elif charge_type in esp_charges:
        name = 'ESP'
    else:
        raise NotImplementedError("Charge type '{0}' is not implemented."
                                  .format(charge_type))
    return ' ' + name + ' charges:'
