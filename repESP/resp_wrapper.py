"""Wrappers for performing ESP fitting with the ``resp`` program"""

from repESP.charges import Charge
from repESP.esp_util import EspData, write_resp_esp
from repESP.equivalence import Equivalence
from repESP.respin_generation import prepare_respin
from repESP.respin_generation import RespStage1RespinGenerator, RespStage2RespinGenerator
from repESP.respin_generation import FitHydrogensOnlyRespinGenerator, FrozenAtomsRespinGenerator
from repESP.respin_generation import EquivalenceOnlyRespinGenerator
from repESP.respin_format import Respin, write_respin
from repESP.resp_charges_format import write_resp_charges, parse_resp_charges
from repESP.types import Atom, Molecule

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
        write_respin(f, respin)

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
    """Run the ``resp`` program with the given "respin" instructions

    Parameters
    ----------
    esp_data : EspData
        Object containing the data normally provided in the .esp file i.e. atom
        coordinates and ESP field values at the points to be used in the fitting.
    respin : Respin
        Instructions for the fitting in the `resp` program input format.
    initial_charges : Optional[typing.List[Charge]], optional
        Initial charges to be used for the fitting, which may be required by
        the fitting method or may simply be provided as an initial guess for
        the fitting algorithm. `ValueError` will be raised if the `respin`
        argument specifies that initial charges are expected (`iqopt` not equal
        to 1) but this argument is not provided. Defaults to None.
    generate_esout : bool, optional
        Whether to produce an "esout" file containing the ESP field values at
        the fitting points, reproduced from the fitted charges. This can also
        be achieved with the `repESP` library and thus shouldn't be necessary.
        Defaults to False. If this option is set to True, the user should also
        set the `save_intermediates_to` option, otherwise the "esout" file will
        not be accessible.
    save_intermediates_to : Optional[str], optional
        Which directory to save the ``resp`` program output to. The directory
        must not exist and will be created by the program. If the user does
        not require access to the files, this option should remain at the default
        of None. The calculation will then be run in the OS temporary directory.

    Returns
    -------
    typing.List[Charge]
        The charges fitted using the ``resp`` program according to the "respin"
        instructions.
    """

    if save_intermediates_to is None:
        with tempfile.TemporaryDirectory() as temp_dir_name:
            return _run_resp_in_dir(esp_data, respin, initial_charges, generate_esout, temp_dir_name)
    else:
        os.mkdir(save_intermediates_to)
        return _run_resp_in_dir(esp_data, respin, initial_charges, generate_esout, save_intermediates_to)


# NOTE: If alternative interface is to be implemented in place of the two
# respin files, try the `variants` library presented by Paul Ganssle.
def run_two_stage_resp(
    esp_data: EspData,
    respin1: Respin,
    respin2: Respin,
    initial_charges: Optional[List[Charge]]=None,
    generate_esout: bool=False,
    save_intermediates_to: Optional[str]=None
) -> List[Charge]:
    """Apply the two-stage procedure to fit RESP charges

    This function takes the two "respin" files, as this is the simplest way to
    provide information about equivalence as well as methyl and methylene groups.
    Another interface is possible where this information is provided directly
    and may be provided in the future.

    Parameters
    ----------
    esp_data : EspData
        See `run_resp` function parameter
    respin1 : Respin
        Instructions for 1st stage RESP fitting.
    respin2 : Respin
        Instructions for 2nd stage RESP fitting.
    initial_charges : Optional[typing.List[Charge]], optional
        See `run_resp` function parameter
    generate_esout : bool, optional
        See `run_resp` function parameter
    save_intermediates_to : Optional[str], optional
        See `run_resp` function parameter

    Returns
    -------
    typing.List[Charge]
        The fitted two-stage RESP charges.
    """
    # Some verification of respin file compatibility should be performed.
    total_charge = respin1.charge
    molecule = respin1.molecule

    get_calc_dir = lambda stage: (
        f"{save_intermediates_to}/stage_{stage}"
        if save_intermediates_to is not None else None
    )

    respin1_generated = prepare_respin(
        RespStage1RespinGenerator(respin1.ivary),
        total_charge,
        molecule,
        read_charges=initial_charges is not None
    )

    resp1_charges = run_resp(
        esp_data,
        respin1,
        initial_charges,
        generate_esout,
        get_calc_dir(1)
    )

    respin2_generated = prepare_respin(
        RespStage2RespinGenerator(respin2.ivary),
        total_charge,
        molecule,
        read_charges=True
    )

    return run_resp(
        esp_data,
        respin2,
        resp1_charges,
        generate_esout,
        get_calc_dir(2)
    )


# NOTE: It may be more natural (i.e. in line with this library's datastructures)
# to use Field and Molecule[AtomWithCoords] instead of EspData. Currently,
# the atom coordinates are hidden in EspData, making it less clear what the
# input of the fitting is.
def fit_with_equivalencing(
    esp_data: EspData,
    equivalence: Equivalence,
    molecule: Molecule[Atom],
    total_charge: int,
    initial_charges: Optional[List[Charge]]=None,
    generate_esout: bool=False,
    save_intermediates_to: Optional[str]=None
) -> List[Charge]:
    """Fit charges to the provided ESP subject to equivalence relations

    Parameters
    ----------
    esp_data : EspData
        See `run_resp` function parameter
    equivalence : Equivalence
        The chemical equivalence relations between atoms of the molecule.
    molecule : Molecule[Atom]
        The molecule which charges are being fitted. Only atom identities are required.
    total_charge : int
        The total charge of the molecule.
    initial_charges : Optional[typing.List[Charge]], optional
        See `run_resp` function parameter
    generate_esout : bool, optional
        See `run_resp` function parameter
    save_intermediates_to : Optional[str], optional
        See `run_resp` function parameter

    Returns
    -------
    typing.List[Charge]
        The charges fitted to the provided ESP subject to equivalence relations.
    """

    respin_generator = EquivalenceOnlyRespinGenerator(equivalence)
    respin = prepare_respin(
        respin_generator,
        total_charge,
        molecule,
        read_charges=initial_charges is not None
    )

    return run_resp(
        esp_data,
        respin,
        initial_charges,
        generate_esout,
        save_intermediates_to
    )


def fit_hydrogens_only(
    esp_data: EspData,
    equivalence: Equivalence,
    molecule: Molecule[Atom],
    total_charge: int,
    initial_charges: List[Charge],
    generate_esout: bool=False,
    save_intermediates_to: Optional[str]=None
) -> List[Charge]:
    """Fit hydrogen atom charges to the provided ESP subject to equivalence relations

    Parameters
    ----------
    esp_data : EspData
        See `run_resp` function parameter
    equivalence : Equivalence
        The chemical equivalence relations between atoms of the molecule.
    molecule : Molecule[Atom]
        The molecule which charges are being fitted. Only atom identities are required.
    total_charge : int
        The total charge of the molecule.
    initial_charges : typing.List[Charge]
        Initial charges to be used for the fitting. Charges on atoms other than
        hydrogen atoms will be fixed at the provided values.
    generate_esout : bool, optional
        See `run_resp` function parameter
    save_intermediates_to : Optional[str], optional
        See `run_resp` function parameter

    Returns
    -------
    typing.List[Charge]
        The hydrogen charges fitted to the provided ESP subject to equivalence
        relations. Other charges have values fixed at values from `initial_charges`.
    """

    respin_generator = FitHydrogensOnlyRespinGenerator(equivalence, molecule)
    respin = prepare_respin(
        respin_generator,
        total_charge,
        molecule,
        read_charges=True
    )

    return run_resp(
        esp_data,
        respin,
        initial_charges,
        generate_esout,
        save_intermediates_to
    )


def fit_with_frozen_atoms(
    esp_data: EspData,
    equivalence: Equivalence,
    molecule: Molecule[Atom],
    frozen_atoms: List[int],
    total_charge: int,
    initial_charges: List[Charge],
    generate_esout: bool=False,
    save_intermediates_to: Optional[str]=None
) -> List[Charge]:
    """Fit hydrogen atom charges to the provided ESP subject to equivalence relations

    Parameters
    ----------
    esp_data : EspData
        See `run_resp` function parameter
    equivalence : Equivalence
        The chemical equivalence relations between atoms of the molecule.
    molecule : Molecule[Atom]
        The molecule which charges are being fitted. Only atom identities are required.
    frozen_atoms : List[int]
        List of zero-based indices of atoms in the molecule which charges should
        not be fitted but fixed at the values provided in `initial_charges` argument.
    total_charge : int
        The total charge of the molecule.
    initial_charges : typing.List[Charge]
        Initial charges to be used for the fitting. Charges on atoms specified
        in the `frozen_atoms` arguments will be fixed at the provided values.
    generate_esout : bool, optional
        See `run_resp` function parameter
    save_intermediates_to : Optional[str], optional
        See `run_resp` function parameter

    Returns
    -------
    typing.List[Charge]
        The charges fitted to the provided ESP subject to equivalence relations,
        except where the charges were frozen at initial values.
    """

    respin_generator = FrozenAtomsRespinGenerator(equivalence, frozen_atoms)
    respin = prepare_respin(
        respin_generator,
        total_charge,
        molecule,
        read_charges=True
    )

    return run_resp(
        esp_data,
        respin,
        initial_charges,
        generate_esout,
        save_intermediates_to
    )
