from repESP import resp, charges
from repESP.field_comparison import rms_and_rep

import shutil
import os
import math
import copy
import pickle


class FitCalc(object):
    """An object creating Gaussian input for modified ESP fitting"""

    def __init__(self, path, basename, method, charge_type, net_charge,
                 multiplicity, iop41=0, iop42=0, iop43=0):
        if path[-1] != '/':
            path += '/'
        self.path, self.basename = path, basename
        self.charge_type, self.method = charge_type, method
        self.iop41, self.iop42, self.iop43 = iop41, iop42, iop43
        self.net_charge, self.multiplicity = net_charge, multiplicity

        self._set_filename()
        self._check_path()

    def _check_path(self):
        if os.path.isfile(self.path):
            raise IOError("The given path is not a direcotry: " + self.path)
        if not os.path.exists(self.path):
            raise IOError("The given path doesn't exist: " + self.path)
        if not os.path.isfile(self.path + self.basename + '.chk'):
            raise IOError("The base checkpoint file {0}.chk doesn't exist in "
                          "directory {1}".format(self.basename + ".chk",
                                                 self.path))
        for fn in ['.com', '.chk', '.log']:
            fn = self.path + self.filename + fn
            if os.path.exists(fn):
                raise FileExistsError("The following file exists: " + fn)

    def _set_filename(self):
        """Set the filename attribute and the IOp string for the `.com` file"""
        result = self.basename + '_' + self.charge_type + '_'
        iop_string = "6/50=1"
        for num in range(41, 44):
            result += "{0:02d}".format(getattr(self, "iop"+str(num)))
            val = getattr(self, "iop" + str(num))
            if val:
                if iop_string:
                    iop_string += ", "
                iop_string += "6/{0}={1}".format(num, val)
        self.filename = result
        self.iop_string = iop_string

    def create_input(self):
        self._check_path()
        shutil.copy(self.path + self.basename + '.chk', self.path +
                    self.filename + '.chk')
        with open(self.path + self.filename + '.com', 'w') as f:
            f.write("%chk={0}.chk\n".format(self.filename))
            f.write("# {0} Geom=Check Pop={1}".format(
                self.method, self.charge_type.upper()))
            if self.iop_string:
                f.write(" IOp({0})".format(self.iop_string))
            f.write("\n\nCreated by repESP script `fit_points`\n\n{0} {1}"
                    .format(self.net_charge, self.multiplicity))
            f.write("\n\n{0}.esp".format(self.filename))


class IOpCalcSet(object):
    """Object creating sets of IOp values from given ranges and options"""

    def __init__(self, iop41=None, iop42=None, iop43=None,
                 thickness=None, thick_round='closest'):
        """Creatre object for creating sets of IOp values

        Parameters `iop41`, `iop42` and `iop43` can be specified as `None`, 0,
        an integer or a list of integers. If `None` or 0 is given, the default
        value will be looked-up for `iop41` and `iop42`; for `iop43` 0 will be
        passed, as it requests an option which cannot be specified otherwise.
        At least one of the parameters must be in the form of a list, as its
        length will be used as the requested length of other lists. If more
        than one of the three parameters is given as a list, they must be of
        equal lengths. If a parameter is given as an integer, a list will be
        created consisting of the same element repeated.

        The default value of `iop43` requests Gaussian's own formula
        calculating it based on the number of layers (`iop41`). However, to
        keep the overall thickness of the set of points, another formula must
        be used. This an be achieved by leaving `iop43` as None and specifying
        the desired thickness instead. To use the default thickness (resulting
        from Gaussian's combination of default number of layers and increment
        formula), specify a thickness of:

        .. math:: (#l-1)*0.4/sqrt(#l) = (4-1)*0.4/sqrt(4) = 0.3

        However, since `iop43` can only be an integer (which is then
        multiplied by 0.01 to obtain the actual value of the increment), it
        will not be possible to keep a constant value of thickness. The default
        method of rounding the increment is 'closest' but you can also select
        'up' or 'down' to obtain a range on your plot.
        """
        if thickness is not None and iop43 is not None:
            raise ValueError("Option `thickness` cannot be set if `iop43` is "
                             "not `None`.")

        self.iop41, self.iop42, self.iop43 = iop41, iop42, iop43
        self.thickness = thickness
        self.thick_round = self._get_rounding_function(thick_round)
        self._check_iop_lists()
        self._check_iop41(iop41)

    def create_param_list(self):
        """Creating sets of IOp values from given ranges and options

        The resulting lists can be iterated together or a meshgrid can be
        formed, depending on the intention.
        """
        result = []
        for iop, default in (('iop41', 4), ('iop42', 1), ('iop43', 0)):
            iop = getattr(self, iop)
            if iop is None or iop == 0:
                iop = default
            if type(iop) is not list:
                iop = [iop]*self.iop_len
            result.append(iop)

        if self.thickness is not None:
            # Calculate intervals according to thickness
            result[-1] = [self.thick_round(100*self.thickness/(layers_count-1))
                          for layers_count in result[0]]

        self._check_result(result)
        return result

    def _check_result(self, result):
        for result_elem in result:
            for elem in result_elem:
                if elem > 99:
                    raise ValueError(
                        "One of the calculated IOp values is greater than 99, "
                        "which is not implemented as it's probably not "
                        "intentional. Try changing input values.")

    @staticmethod
    def _get_rounding_function(thick_round):
        if thick_round == 'closest':
            return round
        elif thick_round == 'down':
            return math.floor
        elif thick_round == 'up':
            return math.ceil
        else:
            raise NotImplementedError("The rounding function '{0}' is not "
                                      "implemented.".format(thick_round))

    def _check_iop_lists(self):
        given_lists = [len(elem) for elem in (
            self.iop41, self.iop42, self.iop43) if type(elem) is list]
        if not len(given_lists):
            raise ValueError("At least one IOp needs to be given as a list.")
        # Check if those which are lists, are of the same lengths
        # http://stackoverflow.com/a/3844948
        if not given_lists.count(given_lists[0]) == len(given_lists):
            raise ValueError("The given IOp lists are of unequal lengths: " +
                             given_lists)
        self.iop_len = given_lists[0]

    @staticmethod
    def _check_iop41(iop41):
        if type(iop41) is not list:
            iop41 = [iop41]
        for elem in iop41:
            if elem is not None and elem < 4:
                raise ValueError("The number of layers needs to be greater or "
                                 "equal to 4.")


charge_type = 'mk'
path = '../data/methane/fit_points-1/'

# PART 1
if True:
    os.mkdir(path)
    shutil.copy('../data/methane/input/methane.chk', path)

    # A simple investigation --- varying only IOp42
    calc_set = IOpCalcSet(iop42=list(range(1, 6)))
    calcs = []
    print(calc_set.create_param_list())
    for iop41, iop42, iop43 in zip(*calc_set.create_param_list()):
        calc = FitCalc(path, 'methane', 'MP2/6-31+G(d,p)', charge_type, 0, 1,
                       iop41, iop42, iop43)
        calc.create_input()
        calcs.append(calc)
        print("Created file: ", calc.filename)

    # Pickle filenames for part 2
    with open("fit_points-calcs.p", 'wb') as f:
        pickle.dump([elem.filename for elem in calcs], f)

# PART 2 --- run when the Gaussian calculations have been completed
if False:
    with open("fit_points-calcs.p", 'rb') as f:
        calcs = pickle.load(f)

    for calc in calcs:
        g = resp.G09_esp(path + calc + '.esp')
        charges.update_with_charges(charge_type, path + calc + '.log',
                                    g.molecule)
        for atom in g.molecule:
            atom.print_with_charge(charge_type)
        min_rms, rep_esp_field = rms_and_rep(g.field, g.molecule, charge_type)
        print(min_rms)
        # More sophisticated analysis of RMS and charges can now follow
