"""Custom exceptions raised by this library"""

class InputFormatError(Exception):
    """Exception raised when given input does not match expected format

    This error will only be raised when parsing external input files, e.g.
    Gaussian output files or ``resp`` input files. Other malformed input from
    the user will result in a `ValueError` being raised, as is idiomatic in
    Python.
    """
    pass
