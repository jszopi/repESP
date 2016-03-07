import shutil
import os
import math


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
        """Create object for creating sets of IOp values

        .. This is not a finished comprehensive docstring. In partiular, the
            meaning of the IOps is not explained. I have only polished the bits
            pertaining to the tricky ``thickness`` calculation.

        Parameters ``iop41``, ``iop42`` and ``iop43`` can be specified as
        ``None``, 0, an integer or a list of integers. If ``None`` or 0 is
        given, the default value will be looked-up for ``iop41`` (4) and
        ``iop42`` (1); for ``iop43`` 0 will be passed, as it requests an option
        which cannot be specified otherwise (Gaussian's own formula, see
        :ref:`below <gaussian_formula>`).

        At least one of the parameters must be in the form of a list,
        as its length will be used as the requested length of other lists. If
        more than one of the three parameters is given as a list, they must be
        of equal lengths. If a parameter is given as an integer, a list will be
        created consisting of the same element repeated.

        The default value of ``iop43`` requests Gaussian's own formula, which
        calculates the shell interval (:math:`d`) based on the number of layers
        (``iop41``, :math:`n_l`):

        .. _gaussian_formula:
        .. math::

            d = \\frac{0.4}{\sqrt{n_l}}

        For the default number of shells (4), the interval is 0.2. Note that,
        according to the MK method, the interval between shells is specified in
        multiples of van der Waals radii.

        .. note:: The original MK paper (Singh1984) cites the MS program by M.
            Connolly for the data on van der Waals (vdW) radii. However, as I
            found when investigating the fitting points around methane, the
            value used by Gaussian for the vdW radius of hydrogen is 1.53 Å,
            which is far off the usual value around 1.2 Å. I suppose that the
            carbon radius used in the Conolly program (`link`_) may have been
            mistakenly used.

            .. _link: http://sites.bmb.uga.edu/wampler/8200/size/MS-radii.html

            .. The parent link is:
                http://sites.bmb.uga.edu/wampler/8200/size/connolly.html#SURF

            Interestingly, my investigation of methane also showed that the
            radius used to exclude the inner points (i.e. calculate the first
            shell at 1.4 vdW radius) is different. It appears to be around 1.18
            Å, so more in line with the usual vdW radius for hydrogen.

        One may also wish to keep the overall thickness of the set of points
        constant. This cannot be achieved through IOps, since the increment
        between MK shells is in multiples of van der Waals radii, so is
        variable depending on the region of the molecule. Therefore thickness
        can be set assuming a single van der Waals radius, that of hydrogen.
        This can be achieved by leaving ``iop43`` as None and specifying the
        desired thickness instead with the ``'thickness'`` keyword. Gaussian's
        Merz-Singh-Kollman radius for hydrogen (1.53 Å) will then be used. The
        thickness resulting from Gaussian's combination of default number of
        layers (4) and the increment calculated with the aforementioned default
        formula is:

        .. math::

           \\newcommand{angstrom}{\\mathrm{\\unicode[serif]{xC5}}}

           1.53 \\angstrom \\times d (n_l-1) = 1.53 \\angstrom \\times (4 - 1)
           \\frac{0.4} {\sqrt{4}} = 0.918 \\angstrom

        Since ``iop43`` can only be an integer (which is then multiplied by
        0.01 to obtain the actual value of the increment), it will not be
        possible to keep exactly a constant value of thickness. The default
        method of rounding the increment is ``'closest'`` but you can also
        select ``'up'`` or ``'down'`` to obtain a range on your plot.

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
            result[-1] = [self.thick_round(100*self.thickness/1.53/(
                          layers_count-1)) for layers_count in result[0]]

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
