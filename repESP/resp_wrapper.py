"""Wrappers for performing ESP fitting with the `resp` program"""

from .charges import Charge
from .esp_util import EspData, write_resp_esp
from .respin_generation import prepare_respin
from .respin_generation import RespStage1RespinGenerator, RespStage2RespinGenerator
from .respin_generation import FitHydrogensOnlyRespinGenerator, FrozenAtomsRespinGenerator
from .respin_generation import EquivalenceOnlyRespinGenerator
from .respin_format import Respin, write_respin, Equivalence
from .resp_charges_format import write_resp_charges, parse_resp_charges
from .types import Atom, Molecule

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

    if save_intermediates_to is None:
        with tempfile.TemporaryDirectory() as temp_dir_name:
            return _run_resp_in_dir(esp_data, respin, initial_charges, generate_esout, temp_dir_name)
    else:
        os.mkdir(save_intermediates_to)
        return _run_resp_in_dir(esp_data, respin, initial_charges, generate_esout, save_intermediates_to)


def run_two_stage_resp(
    esp_data: EspData,
    respin1: Respin,
    respin2: Respin,
    initial_charges: Optional[List[Charge]]=None,
    generate_esout: bool=False,
    save_intermediates_to: Optional[str]=None
) -> List[Charge]:

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


def fit_with_equivalencing(
    esp_data: EspData,
    equivalence: Equivalence,
    molecule: Molecule[Atom],
    total_charge: int,
    initial_charges: Optional[List[Charge]]=None,
    generate_esout: bool=False,
    save_intermediates_to: Optional[str]=None
) -> List[Charge]:

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
