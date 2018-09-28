from repESP.charges import Charge, make_charge
from repESP.esp_util import EspData, parse_gaussian_esp
from repESP.types import *
from repESP.resp_util import run_resp, run_two_stage_resp, fit_hydrogens_only
from repESP.resp_util import fit_with_frozen_atoms, fit_with_equivalencing
from repESP.respin_util import Respin, Equivalence

from my_unittest import TestCase

from io import StringIO


class TestResp(TestCase):

    def setUp(self) -> None:
        with open("data/methane/methane_mk.esp", 'r') as f:
            self.esp_data = EspData.from_gaussian(parse_gaussian_esp(f))

        # First stage RESP
        self.respin = Respin(
            title="File generated for unit tests only (can been removed).",
            cntrl=Respin.Cntrl(
                nmol=1,
                ihfree=1,
                ioutopt=1,
                qwt=0.0005,
            ),
            wtmol=1.0,
            subtitle="Resp charges for organic molecule",
            charge=0,
            iuniq=5,
            atomic_numbers=[6, 1, 1, 1, 1],
            ivary=Respin.Ivary([0, 0, 0, 0, 0])
        )

        self.result_charges = run_resp(
            self.esp_data,
            self.respin
        )

    def test_resp_charges(self) -> None:
        self.assertListsAlmostEqual(
            self.result_charges,
            # Expected values from resp calculations done a while back
            [-0.407205, 0.101907, 0.101695, 0.101695, 0.101907]
        )

    def test_two_stage_resp(self) -> None:

        respin2 = Respin(
            title="File generated for unit tests only (can been removed).",
            cntrl=Respin.Cntrl(
                iqopt=2,
                qwt=0.001,
            ),
            wtmol=1.0,
            subtitle="Resp charges for organic molecule",
            charge=0,
            iuniq=5,
            atomic_numbers=[6, 1, 1, 1, 1],
            ivary=Respin.Ivary([0, 0, 2, 2, 2])
        )

        charges = run_two_stage_resp(
            self.esp_data,
            self.respin,
            respin2
        )

        self.assertListsAlmostEqual(
            charges,
            # Expected values from resp calculations done a while back
            [-0.317454, 0.079364, 0.079364, 0.079364, 0.079364]
        )

    def test_equivalencing(self) -> None:

        charges = fit_with_equivalencing(
            esp_data=self.esp_data,
            equivalence=Equivalence([None, None, 1, 1, 1]),
            atomic_numbers=[6, 1, 1, 1, 1],
            total_charge=0
        )

        self.assertListsAlmostEqual(
            charges,
            # Expected values from resp calculations done a while back
            [-0.500040, 0.125010, 0.125010, 0.125010, 0.125010]
        )

    def test_hydrogen_fitting(self) -> None:

        charges = fit_hydrogens_only(
            esp_data=self.esp_data,
            equivalence=Equivalence([None, None, 1, 1, 1]),
            atomic_numbers=[6, 1, 1, 1, 1],
            total_charge=0,
            initial_charges=[make_charge(x) for x in [-0.5, 0, 0, 0, 0]]
        )

        self.assertListsAlmostEqual(
            charges,
            # Expected values of hydrogen to share 0.5 e equally between them
            # to preserve total charge and net neutral charge.
            [-0.5, 0.125, 0.125, 0.125, 0.125]
        )

    def test_fitting_with_frozen_atoms(self) -> None:

        charges = fit_with_frozen_atoms(
            esp_data=self.esp_data,
            equivalence=Equivalence([None, None, 1, 1, 1]),
            atomic_numbers=[6, 1, 1, 1, 1],
            frozen_atoms=[0, 2, 4],
            total_charge=0,
            initial_charges=[make_charge(x) for x in [-0.5, 0, 0, 0, 0]]
        )

        self.assertListsAlmostEqual(
            charges,
            # Expected values of two hydrogen to share 0.5 e equally between them
            # to preserve total charge and net neutral charge.
            [-0.5, 0.25, 0.0, 0.25, 0.0]
        )
