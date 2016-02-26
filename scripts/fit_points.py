import shutil
import os


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
