from fortranformat import FortranRecordWriter
import textwrap
import os

from .cube_helpers import InputFormatError, Atom, Molecule, Field
from .cube_helpers import angstrom_per_bohr
from . import charges


class DuplicateEntryError(Exception):
    pass


class InputValueError(Exception):
    pass


class G09_esp(object):

    def __init__(self, fn, coords_in_bohr=True, allow_dupes=False):
        self._read_in(fn, coords_in_bohr, allow_dupes)

    def _read_in(self, fn, coords_in_bohr, allow_dupes):
        with open(fn, 'r') as f:
            # Checks two first lines
            self._read_header(fn, f)
            self._read_atoms(f, coords_in_bohr)
            self._read_moments(f)
            self._read_esp_points(f, coords_in_bohr, allow_dupes)

    def _read_header(self, fn, f):
        line = f.readline().rstrip('\n')
        if line != " ESP FILE - ATOMIC UNITS":
            raise InputFormatError(
                "The input file {0} does not seem to be the G09 .esp format. "
                "Generate by specifying Pop=MK/CHelp(G) with IOp(6/50=1)"
                .format(fn))
        line = f.readline().split()
        self.charge = int(line[2])
        self.multip = int(line[-1])

    def _read_atoms(self, f, coords_in_bohr):
        line = f.readline().split()
        atom_count = int(line[-1])
        self.molecule = Molecule(self)
        for i in range(atom_count):
            line = f.readline().split()
            identity = line[0]
            atomic_no = Atom.inv_periodic[identity]
            coords = [float(coord.replace('D', 'E')) for coord in line[1:4]]
            # Neglect the ESP value at atoms, which is given by last value
            self.molecule.append(Atom(i+1, atomic_no, coords, coords_in_bohr))

    def _read_moments(self, f):
        assert f.readline().rstrip('\n') == " DIPOLE MOMENT:"
        # Currently not implemented, the lines are just skipped
        for i in range(4):
            f.readline()

    def _read_esp_points(self, f, coords_in_bohr, allow_dupes):
        line = f.readline().split()
        expected = "ESP VALUES AND GRID POINT COORDINATES. #POINTS ="
        assert ' '.join(line[:-1]) == expected
        expected_points_count = int(line[-1])

        points_coords = []
        values = []
        for line in f:
            line = [val.replace('D', 'E') for val in line.split()]
            points_coords.append(tuple(line[1:4]))
            values.append(float(line[0]))

        if len(points_coords) != expected_points_count:
            raise InputFormatError(
                "The number of ESP points {0} does not agree with that "
                "specified at the top of the input file: {1}".format(
                    len(points_coords), expected_points_count))

        try:
            points = Points(points_coords, coords_in_bohr, allow_dupes)
        except DuplicateEntryError:
            raise InputFormatError(
                "Duplicate points in the input file. This might be an artefact"
                " of the algorithm which produced the points. If these points "
                "are to be counted twice, the NonGridField needs to be called "
                "with `allow_dupes=True`")
        except InputValueError as e:
            # Translate the errors when creating Points to errors due to input
            # file format
            raise InputFormatError(e)

        self.field = NonGridField(values, points, 'esp', field_info=['input'])


class Points(list):

    def __init__(self, points_coords, coords_in_bohr=False, allow_dupes=True):
        super().__init__()
        self.allow_dupes = allow_dupes
        if not self.allow_dupes:
            self.points_dict = {}

        for point_coords in points_coords:
            self.append(self._check_and_create_point(point_coords,
                                                     coords_in_bohr))

    def _check_and_create_point(self, point_coords, coords_in_bohr):
        if len(point_coords) != 3:
            raise InputValueError(
                "Encountered a point with a number of coordinates {0}, which "
                "is different from 3.".format(len(point_coords)))

        if not self.allow_dupes:
            for point_coord in point_coords:
                if type(point_coord) is not str:
                    raise InputValueError(
                        "If no duplicates are allowed, `points` must be"
                        " a list of tuples of *strings*. Encountered type {0} "
                        "instead.".format(type(point_coord)))

            if point_coords in self.points_dict:
                raise DuplicateEntryError("Encountered a duplicate point: {0}"
                                          .format(point_coords))

        try:
            result = [float(point_coord) for point_coord in point_coords]
        except ValueError:
            raise InputValueError("Couldn't convert coordinates {0} to float."
                                  .format(point_coords))

        if coords_in_bohr:
            result = [angstrom_per_bohr*point_coord for point_coord in result]

        return result


class NonGridField(Field):

    def __init__(self, values, points, field_type, field_info=None):
        """Create a NonGridField from given coordinates and values

        Parameters
        ----------
        points : Points

        values : List
            The list of values at *corresponding* coordinates.
        """
        super().__init__(values, field_type, field_info, check_nans=False)
        self.points = points

        if type(points) is not Points:
            raise TypeError("Expected type Points for the points argument, "
                            "instead got {0}".format(type(points)))

        if len(points) != len(values):
            raise ValueError("The number of points {0} is different from the "
                             "number of values {1}.".format(len(points),
                                                            len(values)))

    def write_to_file(self, output_fn, molecule, write_coords_in_bohr=True):
        # Numeric formats specified in resp input specification
        # http://upjv.q4md-forcefieldtools.org/RED/resp/#other3
        header_format = FortranRecordWriter('2I5')
        atoms_format = FortranRecordWriter('17X,3E16.7')
        esp_points_format = FortranRecordWriter('1X,4E16.7')

        with open(output_fn, 'x') as f:
            f.write(header_format.write([len(molecule),
                                         len(self.points)]) + "\n")
            for atom in molecule:
                if write_coords_in_bohr:
                    coords = [atom_coord/angstrom_per_bohr for atom_coord in
                              atom.coords]
                else:
                    coords = atom.coords
                f.write(atoms_format.write(coords) + "\n")
            for point_coords, esp_val in zip(self.points, self.values):
                if write_coords_in_bohr:
                    point_coords = [point_coord/angstrom_per_bohr for
                                    point_coord in point_coords]
                f.write(esp_points_format.write([esp_val] + point_coords)+"\n")


def _read_respin(fn, ref_molecule=None):
    with open(fn, 'r') as inp:
        # Rewind to end of cntrl section
        line = inp.readline()
        while "&end" not in line:
            line = inp.readline()
        # Skip two lines ...
        for i in range(3):
            line = inp.readline()
        # ... and the third one will be `charge, iuniq`
        charge, iuniq = [int(elem) for elem in line.split()]

        # Create a molecule for consistency checks. It will not be returned, so
        # should be garbage-collected.
        molecule = Molecule(None)
        ivary_list = []
        for i, line in enumerate(inp):
            if len(line.split()) != 2:
                break
            atom = Atom(i+1, int(line.split()[0]))
            # Crucial bit: reading in ivary
            ivary_list.append(int(line.split()[1]))
            molecule.append(atom)

    # Check input file consistency
    if len(molecule) != iuniq:
        raise InputFormatError("The number of atoms {0} doesn't agree with"
                               " the `iuniq` value in the input file: {1}"
                               .format(len(molecule), iuniq))
    # Check the molecule against a reference molecule
    if ref_molecule is not None and molecule != ref_molecule:
        molecule.verbose_compare(ref_molecule)
        raise InputFormatError("The molecule in the .respin file differs "
                               "from the other molecule as shown above.")

    return ivary_list, charge, iuniq


common_respin_head = """
                     &cntrl

                     nmol = 1,
                     ihfree = 1,
                     ioutopt = 1,
                     """

common_respin_tail = """
                      &end
                         1.0
                     Resp charges for organic molecule
                     """


def _get_respin_content(respin_type, read_input_charges):
    """Check if respin type is implemented and return the input content"""
    if respin_type not in ['1', '2', 1, 2, 'h', 'u']:
        raise ValueError("`respin_type` {0} is not implemented.".format(
                         respin_type))

    result = "RESP input of type '{0}' generated by the repESP program".format(
        respin_type) + "\n" + textwrap.indent(textwrap.dedent(
            common_respin_head), ' ')

    if str(respin_type) == '1':
        result += " qwt = 0.00050,\n"
    elif str(respin_type) == '2':
        result += " qwt = 0.00100,\n"
    elif str(respin_type) == 'u':
        result += " qwt = 0.00000,\n"

    if read_input_charges:
        result += " iqopt = 2,\n"

    return result + textwrap.dedent(common_respin_tail)


def _write_modified_respin(respin_type, molecule, ivary_list, charge, iuniq,
                           fn_out, read_input_charges):
    with open(fn_out, 'w') as out:
        out.write(_get_respin_content(respin_type, read_input_charges))
        numbers = FortranRecordWriter('2I5')
        # `charge, iuniq` line
        print(numbers.write([charge, iuniq]), file=out)
        for atom, ivary in zip(molecule, ivary_list):
            print(numbers.write([atom.atomic_no, ivary]), file=out)

        print(file=out)


def _check_ivary(check_ivary, molecule, ivary_list):
    if not check_ivary:
        return
    print("\nPlease check if the following RESP input is what you want:\n")
    for atom, ivary in zip(molecule, ivary_list):
        print(atom, end='')
        if ivary == -1:
            print(", frozen")
        elif ivary > 0:
            print(", equivalenced to atom", molecule[ivary-1].label)
        else:
            print()


def _get_input_files(input_dir, respin1_fn, respin2_fn, esp_fn):
    if input_dir[-1] != '/':
        input_dir += '/'
    input_dir_contents = os.listdir(input_dir)
    extensions = ['.respin1', '.respin2', '.esp']
    result = []
    filenames = [respin1_fn, respin2_fn, esp_fn]
    for extension, fn in zip(extensions, filenames):
        candidates = [f for f in input_dir_contents if f.endswith(extension)]
        if not len(candidates):
            raise FileNotFoundError(
                "The input directory {0} doesn't contain any {1} files."
                .format(input_dir, extension))
        if fn in candidates:
            result.append(input_dir + fn)
            continue
        elif len(candidates) > 1:
            raise InputFormatError(
                "{0} {1} files found in the input directory {2}. Please "
                "specify the filename of the file to be used."
                .format(len(candidates), extension, input_dir))
        result.append(input_dir + candidates[0])
    return result


def run_resp(input_dir, calc_dir_path, resp_type='two_stage', inp_charges=None,
             check_ivary=True, respin1_fn=None, respin2_fn=None, esp_fn=None):
    """Runs RESP fitting and returns a molecule updated with resulting charges

    The necessary input files (``.esp``, ``.respin1`` and ``.respin2``) will be
    found in the input directory by extension. A new directory
    ``calc_dir_path`` will be created to keep all the intermediate and output
    files of the ``resp`` program.

    Parameters
    ----------
    input_dir : str
        Directory containing the input files.
    calc_dir_path : str
        Path to the new directory to be created.
    resp_type : {'two_stage', 'h_only', 'unrest'}, optional
        The default ``two_stage`` option requests the normal two-stage RESP
        fitting. The ``ivary`` options are taken unaltered from the two
        original ``.respin`` files. ``unrest`` requests unconstrained
        optimization through one-stage RESP with zero restraint weight
        (``qwt``). ``h_only`` also performs one-stage fitting but freezes all
        atoms except for hydrogens at input values. Both ``h_only`` and
        ``unrest`` read atom equivalence from the ``.respin`` files (``ivary``
        values). To verify that the equivalence is as expected, leave the
        ``check_ivary`` option enabled.

        .. The RESP logic is not explained very well in the papers but I've
        additionally re-engineered the ``resp`` program's logic to be sure that
        reading both the ``respin`` files will give the desired behaviour. In
        fact, it's pretty simple. In the first stage atoms of the methyl and
        methylene groups are free, while all the others are equivalenced. In
        the second stage the former are equivalenced, while all the others are
        frozen.

    inp_charges : List[float], optional
        The input charges. Defaults to ``None``, which causes no ``iqopt``
        command being specified in the ``&cntrl`` section of the ``.respin``
        file. This causes it to default to 1 and 'reset all initial charges to
        zero'.
    check_ivary : bool, optional
        Verbosely report the RESP ``ivary`` actions to be performed by the
        ``resp`` program and ask the user if this is desired behaviour (no
        prompt). This is recommended as the equivalence information is taken
        from the ``.respin`` files, which are generated externally (likely by
        the ``respgen`` program) and, for one-stage RESP (``h_only`` or
        ``unrest``), additionaly modified by this program.
    respin1_fn,respin2_fn,esp_fn : str, optional
        The filenames of input files. These should be specified if there are
        more files with the same extension in the input directory.

    Returns
    -------
    Molecule
        Molecule created based on the ``.esp`` input file, updated with RESP
        charges.

    """
    if calc_dir_path[-1] != '/':
        calc_dir_path += '/'
    os.mkdir(calc_dir_path)
    respin1_fn, respin2_fn, esp_fn = _get_input_files(input_dir, respin1_fn,
                                                      respin2_fn, esp_fn)
    g09_esp = G09_esp(esp_fn)
    molecule = g09_esp.molecule
    if inp_charges is not None and len(inp_charges) != len(molecule):
        raise InputFormatError("The list of input charges is of length {0} but"
                               " the molecule considered has {1} atoms."
                               .format(len(inp_charges), len(molecule)))
    # Create the corrected .esp file
    g09_esp.field.write_to_file(calc_dir_path + "corrected.esp", molecule)
    # Dump the input charges
    if inp_charges is not None:
        charges._update_molecule_with_charges(molecule, inp_charges,
                                              'resp_inp')
        charges.dump_charges_to_qout(molecule, 'resp_inp', calc_dir_path +
                                     "input.qout")

    if resp_type == 'two_stage':
        charges_out_fn = _resp_two_stage(
            calc_dir_path, respin1_fn, respin2_fn, molecule, check_ivary,
            inp_charges is not None)
    elif resp_type == 'h_only' or resp_type == 'unrest':
        charges_out_fn = _resp_one_stage(resp_type[0], calc_dir_path,
                                         respin1_fn, respin2_fn, molecule,
                                         check_ivary, inp_charges is not None)
    else:
        raise ValueError("RESP fitting type '{0}' was not recognized."
                         .format(resp_type))

    # Update the molecule with new 'resp' charges
    print("\nNOTE: The .qout file was generated by the program, so you can "
          "ignore the following warning.")
    charges.update_with_charges('resp', calc_dir_path + charges_out_fn,
                                molecule)

    return molecule


def _resp_one_stage(resp_type, calc_dir_path, respin1_fn, respin2_fn, molecule,
                    check_ivary, read_input_charges):
    """A common function for one-stage RESP ('h_only' and 'unrest')

    Atom equivalence will be taken from the ``.respin`` files (``ivary``
    values).
    """
    # MODIFY THE .RESPIN FILE
    # Read in .respin1
    ivary_list1, charge1, iuniq1 = _read_respin(respin1_fn,
                                                ref_molecule=molecule)
    # Read in .respin2
    ivary_list2, charge2, iuniq2 = _read_respin(respin2_fn,
                                                ref_molecule=molecule)
    assert charge1 == charge2
    assert iuniq1 == iuniq2
    # Modify ivary list and write to new input file
    ivary_list = _modify_ivary_list(resp_type, molecule, ivary_list1,
                                    ivary_list2)
    _check_ivary(check_ivary, molecule, ivary_list)
    _write_modified_respin(resp_type, molecule, ivary_list, charge1, iuniq1,
                           calc_dir_path + "input.respin",
                           read_input_charges=read_input_charges)

    # RUN RESP
    input_charges_option = "-q input.qout " if read_input_charges else ""
    os.system("cd {0}; resp -i input.respin -o output.respout -e "
              "corrected.esp ".format(calc_dir_path) + input_charges_option +
              "-t charges.qout")
    return "charges.qout"


def _modify_ivary_list(resp_type, molecule, ivary_list1, ivary_list2):
    result = []
    for atom, ivary1, ivary2 in zip(molecule, ivary_list1, ivary_list2):
        if resp_type in ['h', 'u']:
            # Extract equivalence from the two default RESP inputs. This simple
            # condition ensures that equivalence (positive number) is picked
            # over free fitting (zero), which is picked over freezing charges
            # (negative number).
            ivary = max(ivary1, ivary2)
        else:
            raise NotImplementedError("Modification of ``ivary`` values not "
                                      "implemented for resp_type '{0}'."
                                      .format(resp_type))
        if resp_type == 'h':
            # Additionally freeze non-hydrogens
            ivary = ivary if atom.atomic_no == 1 else -1
        result.append(ivary)
    return result


def _resp_two_stage(calc_dir_path, respin1_fn, respin2_fn, molecule,
                    check_ivary, read_input_charges):
    """Run the two-stage RESP but with potentially non-zero initial charges"""
    if check_ivary:
        print("\nTwo-stage RESP --- `ivary` values were not modified from the "
              "input `respin` files but you may want to inspect them "
              "nevertheless.\n\nSTAGE 1")

    ivary_list1, charge1, iuniq1 = _read_respin(respin1_fn,
                                                ref_molecule=molecule)
    # ivary_list1 used without modification. Modify the .respin1 file only to
    # ask to load initial charges if `read_input_charges` is True.
    _check_ivary(check_ivary, molecule, ivary_list1)
    _write_modified_respin('1', molecule, ivary_list1, charge1, iuniq1,
                           calc_dir_path + "input1.respin",
                           read_input_charges=read_input_charges)

    if check_ivary:
        print("\nSTAGE 2")
    # Although copying `.respin2` would suffice, _write_modified_respin is
    # called here as well for consistency. The respin file content could
    # potentially differ from that produced by `respgen` if its defaults change
    ivary_list2, charge2, iuniq2 = _read_respin(respin2_fn,
                                                ref_molecule=molecule)
    _check_ivary(check_ivary, molecule, ivary_list2)
    _write_modified_respin('2', molecule, ivary_list2, charge2, iuniq2,
                           calc_dir_path + "input2.respin",
                           read_input_charges=True)

    assert charge1 == charge2
    assert iuniq1 == iuniq2

    # Run resp
    input_charges_option = "-q input.qout " if read_input_charges else ""
    os.system("cd {0}; resp -i input1.respin -o output1.respout -e "
              "corrected.esp ".format(calc_dir_path) + input_charges_option +
              "-t charges1.qout -s esout1 -p punch1")
    os.system("cd {0}; resp -i input2.respin -o output2.respout -e "
              "corrected.esp -q charges1.qout -t charges2.qout -s esout2 "
              "-p punch2".format(calc_dir_path))
    return "charges2.qout"
