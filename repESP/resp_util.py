from .charges import Charge, make_charge
from .esp_util import EspData, write_resp_esp
from .respin_util import Respin, _write_respin
from .resp_charges_util import write_resp_charges, parse_resp_charges

import os
import subprocess
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
