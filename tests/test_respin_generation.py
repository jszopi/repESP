from repESP.types import *
from repESP.respin_util import Respin, Equivalence
from repESP.respin_generation import prepare_respin
from repESP.respin_generation import RespStage1RespinGenerator, FitHydrogensOnlyRespinGenerator
from repESP.respin_generation import EquivalenceOnlyRespinGenerator, FrozenAtomsRespinGenerator

from my_unittest import TestCase

from copy import deepcopy
from io import StringIO


class TestRespinGenerationOnStage1Resp(TestCase):

    def setUp(self) -> None:
        # First stage RESP
        self.respin = Respin(
            title="Respin file prepared by `repESP` to perform Resp Stage 1.",
            cntrl=Respin.Cntrl(
                qwt=0.0005
            ),
            wtmol=1.0,
            subtitle="Resp charges for organic molecule",
            charge=0,
            iuniq=5,
            atomic_numbers=[6, 1, 1, 1, 1],
            ivary=Respin.Ivary([0, 0, 0, 0, 0])
        )

    def test_defaults(self) -> None:

        generated = prepare_respin(
            respin_generator=RespStage1RespinGenerator(self.respin),
            total_charge=0,
            atomic_numbers=self.respin.atomic_numbers
        )

        self.assertAlmostEqualRecursive(self.respin, generated)

    def test_with_read_charges(self) -> None:

        generated = prepare_respin(
            respin_generator=RespStage1RespinGenerator(self.respin),
            total_charge=0,
            atomic_numbers=self.respin.atomic_numbers,
            read_charges=True
        )

        expected = deepcopy(self.respin)
        expected.cntrl.iqopt = 2

        self.assertAlmostEqualRecursive(expected, generated)

    def test_with_total_charge(self) -> None:

        generated = prepare_respin(
            respin_generator=RespStage1RespinGenerator(self.respin),
            total_charge=1,
            atomic_numbers=self.respin.atomic_numbers
        )

        expected = deepcopy(self.respin)
        expected.charge = 1

        self.assertAlmostEqualRecursive(expected, generated)

    def test_with_title_and_subtitles(self) -> None:

        generated = prepare_respin(
            respin_generator=RespStage1RespinGenerator(self.respin),
            total_charge=0,
            atomic_numbers=self.respin.atomic_numbers,
            title="Custom title",
            subtitle="Custom subtitle"
        )

        expected = deepcopy(self.respin)
        expected.title = "Custom title"
        expected.subtitle = "Custom subtitle"

        self.assertAlmostEqualRecursive(expected, generated)


class TestFitingHydrogensOnly(TestCase):

    def test_fitting_hydrogens_only(self) -> None:

        respin = Respin(
            title="Respin file prepared by `repESP` to perform fitting of hydrogen atoms.",
            cntrl=Respin.Cntrl(
                iqopt=2,
                ihfree=0,
                qwt=0.0
            ),
            wtmol=1.0,
            subtitle="Resp charges for organic molecule",
            charge=0,
            iuniq=5,
            atomic_numbers=[6, 1, 1, 1, 1],
            ivary=Respin.Ivary([-1, 0, 2, 2, 2])
        )

        generated = prepare_respin(
            respin_generator=FitHydrogensOnlyRespinGenerator(
                Equivalence([None, None, 1, 1, 1]),
                respin.atomic_numbers
            ),
            total_charge=0,
            atomic_numbers=respin.atomic_numbers
        )

        self.assertAlmostEqualRecursive(respin, generated)


class TestEquivalencing(TestCase):

    def test_equivalencing(self) -> None:

        respin = Respin(
            title="Respin file prepared by `repESP` to perform atom equivalencing.",
            cntrl=Respin.Cntrl(
                ihfree=0,
                qwt=0.0
            ),
            wtmol=1.0,
            subtitle="Resp charges for organic molecule",
            charge=0,
            iuniq=5,
            atomic_numbers=[6, 1, 1, 1, 1],
            ivary=Respin.Ivary([0, 0, 2, 2, 2])
        )

        generated = prepare_respin(
            respin_generator=EquivalenceOnlyRespinGenerator(
                Equivalence([None, None, 1, 1, 1])
            ),
            total_charge=0,
            atomic_numbers=respin.atomic_numbers
        )

        self.assertAlmostEqualRecursive(respin, generated)


class TestFittingWithFrozenAtoms(TestCase):

    def test_fitting_with_frozen_atoms(self) -> None:

        respin = Respin(
            title="Respin file prepared by `repESP` to perform fitting with selected atom charges frozen.",
            cntrl=Respin.Cntrl(
                iqopt=2,
                ihfree=0,
                qwt=0.0
            ),
            wtmol=1.0,
            subtitle="Resp charges for organic molecule",
            charge=0,
            iuniq=5,
            atomic_numbers=[6, 1, 1, 1, 1],
            ivary=Respin.Ivary([-1, 0, 2, -1, 2])
        )

        generated = prepare_respin(
            respin_generator=FrozenAtomsRespinGenerator(
                Equivalence([None, None, 1, 1, 1]),
                frozen_atoms=[0, 3]
            ),
            total_charge=0,
            atomic_numbers=respin.atomic_numbers
        )

        self.assertAlmostEqualRecursive(respin, generated)
