from .charges import Charge, make_charge
from .esp_util import EspData, write_resp_esp
from .respin_util import Respin, _write_respin
from .resp_charges_util import write_resp_charges, parse_resp_charges
from .util import get_symbol

from dataclasses import dataclass
from itertools import zip_longest
import os
import subprocess
import sys
from typing import List, Optional
import tempfile


def _run_resp_in_dir(
    esp_data: EspData,
    respin: Respin,
    initial_charges: Optional[List[Charge]],
    generate_esout: bool,
    calc_dir: str,
) -> List[Charge]:

    if respin.cntrl.iqopt in [2, 3] and initial_charges is None:
        raise ValueError("`resp` expected initial charges (`iqopt` is not 1) but none given.")

    respin_fn = "input.respin"
    qin_fn = "charges.qin"
    qout_fn = "charges.qout"
    espot_fn = "espot.esp"

    get_full_path = lambda fn: f"{calc_dir}/{fn}"

    with open(get_full_path(respin_fn), "w") as f:
        _write_respin(f, respin)

    if initial_charges is not None:
        with open(get_full_path(qin_fn), "w") as f:
            write_resp_charges(f, initial_charges)

    with open(get_full_path(espot_fn), "w") as f:
        write_resp_esp(f, esp_data)

    try:
        process = subprocess.run(
            [
                "resp",
                "-i", "input.respin",
                "-o", "output.respout",
                *(["-q", "charges.qin"] if initial_charges is not None else []),
                "-t", "charges.qout",
                "-e", "espot.esp",
                *(["-s", "esout.esp"] if generate_esout else []),
            ],
            cwd=calc_dir,
            check=True,  # raises CalledProcessError
            text=True,
            capture_output=True
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"Running `resp` failed with return code {e.returncode} and the following error:\n{e.stderr}."
        )

    with open(get_full_path("stdout.out"), 'w') as f:
        f.write(process.stdout)

    with open(get_full_path("charges.qout")) as f:
        return parse_resp_charges(f)


def run_resp(
    esp_data: EspData,
    respin: Respin,
    initial_charges: Optional[List[Charge]]=None,
    generate_esout: bool=False,
    save_intermediates_to: Optional[str]=None
) -> List[Charge]:

    if save_intermediates_to is None:
        with tempfile.TemporaryDirectory() as temp_dir_name:
            return _run_resp_in_dir(esp_data, respin, initial_charges, generate_esout, temp_dir_name)
    else:
        os.mkdir(save_intermediates_to)
        return _run_resp_in_dir(esp_data, respin, initial_charges, generate_esout, save_intermediates_to)


@dataclass
class Equivalence:
    # Zero-indexed, None if not equivalenced to any atom. Note similarities in
    # implementation with Respin.Ivary - perhaps should be refactored.
    values: List[Optional[int]]

    def __post_init__(self):
        for i, elem in enumerate(self.values):
            if elem is not None and (elem < 0 or elem >= len(self.values)):
                raise ValueError(
                    f"Value number {i} is not valid as equivalence information "
                    f"in a molecule of {len(self.values)}."
                )

    def describe(self, atomic_numbers: Optional[List[int]]=None, file=sys.stdout):
        """Verbosely report the equivalence information."""
        if atomic_numbers is not None and len(atomic_numbers) != len(self.values):
            raise ValueError(
                f"The number of atoms ({len(atomic_numbers)} is not the same "
                f"as the number of equivalence values ({len(self.values)}."
            )

        zipped = zip_longest(self.values, atomic_numbers if atomic_numbers is not None else [])

        for i, (equivalence, atomic_number) in enumerate(zipped):
            identity = get_symbol(atomic_number) if atomic_numbers is not None else None
            id_str = f" ({identity})" if identity is not None else ""
            equivalence_str = f", equivalenced to atom {equivalence+1}" if equivalence is not None else ""
            print(f"Atom{id_str} number {i+1}{equivalence_str}", file=file)

    @classmethod
    def from_ivary(cls, ivary: Respin.Ivary):
        """Get atom equivalence information from an Respin.Ivary object.

        Note: Ivary objects are specific to `resp` program input and thus may
        not provide information about atom equivalence. The `respin` file may
        have been generated to perform any custom fitting by RESP. However,
        `respin` files containing equivalence information may be created using
        the following respgen command: (TODO).
        """
        return cls([
            # Whether `ivary - 1` fulfills other preconditions will be checked in __post_init__
            None if ivary == 0 else ivary - 1
            for ivary in ivary.values
        ])
