import ipdb
from cube_helpers import Atom, Cube


class NotImplementedError(Exception):
    pass


class InputFortmatError(Exception):
    pass


def mulliken(filename, atoms):
    """Updates the atom lists with Mulliken charges.

    The charges should be given by a Gaussian output file (.log or
    .out). In the future, checkpoint and formatted checkpoint
    formats may be supported.
    """
    if filename[-4:] in ['.log', '.out']:
        _mulliken_from_log(filename, atoms)
    else:
        raise NotImplementedError('File extension {0} not supported.'
                                  .format(filename[-4]))


def _mulliken_from_log(filename, atoms):
    """Updates the atom lists with Mulliken charges

    Note that if the output file contains information about Mulliken
    charges in more than one place, only the last one will be used.
    Also, the atom list is assumed to be in order.
    """
    with open(filename, 'r') as f:
        _goto_last_mulliken_in_log(f)
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
        if next_line[:24] != ' Sum of Mulliken charges':
            raise InputFortmatError('Expected end of Mulliken charges ( \'Sum '
                                    'of Mulliken charges\'), instead got: ' +
                                    next_line)

        for atom, charge in zip(atoms, charges):
            atom.charges['mulliken'] = float(charge)


def _goto_last_mulliken_in_log(file_object):
    """Go to the last Mulliken charges output in a log file.

    Based on:
    http://stackoverflow.com/a/620492
    """
    offset = 0
    start = 0

    for line in file_object:
        offset += len(line)
        if line.rstrip('\n') == ' Mulliken charges:':
            start = offset

    file_object.seek(start)
    # Skip the unnecessary line
    file_object.readline()
