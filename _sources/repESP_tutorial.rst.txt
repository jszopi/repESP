Tutorial
========

Say that you came up with a new method of assigning partial charges.
In this method, the charges on heavy atoms are identical to NPA charges and charges on hydrogen atoms are optimized to best fit the molecular ESP field.
This example will walk you through the steps to compute the charges for a trimethylamine cation (NMe\ :sub:`3`\ H\ :sup:`+`\ ) and evaluate fit quality with respect to the actual ESP field.

.. note::
    The full script is available alongside its expected output in the `example` directory.
    All necessary input data files are available in the `data/NMe3H_plus <https://github.com/jszopi/repESP/tree/master/data/NMe3H_plus>`_ directory of the GitHub repo.
    Gaussian input files are currently not available but the analogous input files for a methane molecule are available in the `data/methane/prep <https://github.com/jszopi/repESP/tree/master/data/methane/prep>`_ directory.

1. Say that you want to use the CHelpG fitting points as the sampling mesh for the molecular ESP field.
   Gaussian will generate such file when passed the ``Pop=CHelpG`` and ``IOp(6/50=1)`` options.
   The name of the .esp file must be given as the last line of the input file.
   In this example, the filename was "NMe3H_plus_chelpg.esp".
   The code to parse the file is as follows::

        from repESP.esp_util import parse_gaussian_esp

        with open("NMe3H_plus_chelpg.esp") as f:
            gaussian_esp_data = parse_gaussian_esp(f)

        molecule = gaussian_esp_data.molecule
        esp_field = gaussian_esp_data.field

2. Run a Gaussian calculation with the ``Pop=NPA`` option.
   The NPA charges can then be parsed from the output file as follows::

        from repESP.gaussian_format import get_charges_from_log, NpaChargeSectionParser

        with open("NMe3H_plus_nbo.log") as f:
            npa_charges = get_charges_from_log(f, NpaChargeSectionParser(), molecule)

3. Prepare data regarding chemical equivalence relationships between atoms.
   An Antechamber file (.ac) is required to generate two "respin" files using the ``respgen`` program (available from the Antechamber program suite).
   To generate the files, run in your shell:

   .. code-block:: shell

        antechamber -i NMe3H_plus.log -fi gout -o NMe3H_plus.ac -fo ac
        respgen -i NMe3H_plus.ac -o NMe3H_plus.respin1 -f resp1
        respgen -i NMe3H_plus.ac -o NMe3H_plus.respin2 -f resp2

   The `Equivalence` object can then be generated::

        from repESP.respin_format import parse_respin, get_equivalence_from_two_stage_resp_ivary

        with open("NMe3H_plus_mk.respin1") as respin1_file:
            respin1 = parse_respin(respin1_file)

        with open("NMe3H_plus_mk.respin2") as respin2_file:
            respin2 = parse_respin(respin2_file)

        equivalence = get_equivalence_from_two_stage_resp_ivary(respin1.ivary, respin2.ivary)

   Alternatively, if you're currently unable to generate the "respin" files, `equivalence` can be specified manually::

        from repESP.equivalence import Equivalence

        equivalence = Equivalence([
            None,  # C0
            None,  # H1
            1,     # H2, equivalenced to H1
            1,     # H3, equivalenced to H1
            0,     # C4, equivalenced to C0
            1,     # H5, equivalenced to H1
            1,     # H6, equivalenced to H1
            1,     # H7, equivalenced to H1
            0,     # C8, equivalenced to C0
            1,     # H9, equivalenced to H1
            1,     # H10, equivalenced to H1
            1,     # H11, equivalenced to H1
            None,  # N12
            None   # H13
        ])

4. Perform fitting partial charges on hydrogen atoms.
   This invokes the ``resp`` program, so it has to be in your ``PATH`` variable::

        from repESP.esp_util import EspData
        from repESP.resp_wrapper import fit_hydrogens_only

        esp_data = EspData.from_gaussian(gaussian_esp_data)

        my_charges = fit_hydrogens_only(
            esp_data,
            equivalence,
            molecule,
            total_charge=1,  # cation, +1 total charge
            initial_charges=npa_charges  # will be used for non-hydrogen atoms
        )

5. To reproduce the ESP from your charges, you have to first create a `Molecule` object identical to the one in the `molecule` variable but with your new charges instead of CHelpG::

        from repESP.charges import AtomWithCoordsAndCharge
        from repESP.types import Molecule

        molecule_with_my_charges = Molecule([
            AtomWithCoordsAndCharge(
                atom.atomic_number,
                atom.coords,
                my_charge
            ) for atom, my_charge in zip(molecule.atoms, my_charges)
        ])

   You can then reproduce the ESP field from your partial charges and calculate fit quality::

        from repESP.calc_fields import esp_from_charges, calc_rms_error, calc_relative_rms_error

        reproduced_esp = esp_from_charges(esp_field.mesh, molecule_with_my_charges)
        rms = calc_rms_error(esp_field.values, reproduced_esp.values)
        rrms = calc_relative_rms_error(esp_field.values, reproduced_esp.values)

   You could also evaluate the fit quality of the original CHelpG and NPA charges for comparison.
   Extending the script to do that is left as an exercise to the reader.
   You should find that the fit quality of your charges is better than that of NPA but poorer than that of CHelpG.
