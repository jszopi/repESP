[![License: GPL v3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.org/jszopi/repESP.svg?branch=master)](https://travis-ci.org/jszopi/repESP)

# What is `repESP`?
`repESP` is a program used in computational chemistry to investigate different methods of assigning partial charges.
It takes the name from reproducing the (molecular) electrostatic potential (ESP) from partial charges.
Partial charges can be extracted from the output of Gaussian 09 or provided by the user.
Current capabilities of the program include:

* comparing the reproduced ESP to the 'true' ESP (operations on `cube` files)
* evaluating fit quality on a provided set of fitting points
* visualization of fitting points around a molecule, coloured by ESP value or error
* averaging (or, for ESP-based methods, equivalencing) charges of symmetry-related and rapidly-exchanging atoms

The program was developed as part of Jan Szopinski's MSci project at the [Hunt Research Group](http://www.huntresearchgroup.org.uk/), Imperial College London.
The program was used to investigate and develop the following novel concepts:

* charge 'flexibility'
* 'adjusted rational' charges

# How to use `repESP`

`repESP` implemented a lot of functionalities, which can be used to examine the various methods of assigning partial charges.
Many of these functionalities are exposed to users as command line (CLI) scripts and are easy to use.
The scripts perform well-defined, modular tasks and can be woven together to create one's own work flow (Unix style).
The following chart presents some of the `repESP` capabilities:

![Data flow](https://github.com/jszopi/repESP/blob/master/docs/diagrams/data_flow.png)

The monospaced text above arrows shows the names of the relevant command line scripts.

When the out-of-the-box functionalities are not sufficient or flexible enough, users can create their own scripts by using the repESP library.

## Documentation

### Scripts

The overview of the available scripts can be found [here](scripts/README.md) and the command line interface help can be found [here](scripts/detailed.md).

### Library

> **Note:** Work on the library is currently ongoing and documentation will be filled in over the coming week or two, leading to a release of version 0.2.0.
> Missing functionalities will be prioritized and added in subsequent revisions.

Documentation for the library isn't currently hosted and must be generated locally by the user:
After installing `repESP`, run the following commands from the `docs` directory:

```sh
sphinx-apidoc -fMe -o source ../repESP
make html
```

The documentation can then be found in the `build/html` subdirectory.
You can open the `index.html` or `modules.html` file in any browser and navigate from there.

# Run `repESP`

The following instructions are for Unix-like systems, like MacOS and Linux.
There shouldn't be anything stopping `repESP` from running on Windows but that hasn't been tested.

## Prerequisites

`repESP` needs input from computational chemistry software to obtain the ESP and partial charges.
Currently only the output from Gaussian 09 is known to be recognized, please open issues for other Gaussian versions or other software.
Further, to run `repESP` you need Python 3.7 and two other free programs, which need to be in your `PATH`.

### `resp` 

`resp` is a free program necessary to perform ESP-fitting.
Download link can be found at the bottom of this [documentation page](http://upjv.q4md-forcefieldtools.org/RED/resp/).
To install you'll have to compile the Fortran source and put the resulting binary in your `PATH`.

### AmberTools

`resp` and `repESP` need information about atom equivalence in the form of `.respin` files.
AmberTools is a suite of free programs ([link](http://ambermd.org/AmberTools16-get.html)).
However, it requires compilation from C, C++ and Fortran sources, which proved tricky on MacOS.

To create the necessary input for a methane molecule, first create the `.ac` file from Gaussian `.log` using `antechamber`:

`antechamber -i methane.log -fi gout -o methane.ac -fo ac`

Then run `respgen` twice:

```sh
respgen -i methane.ac -o methane.respin1 -f resp1
respgen -i methane.ac -o methane.respin2 -f resp2
```

## Install

Installing `repESP` is standard and simple:

1. Choose and go to an installation directory (e.g. `~/bin`)
2. Get the code: `git clone https://github.com/jszopi/repESP`
3. From inside the repESP directory run: `pip3 install .`

With the following caveats regarding the final point:

* While version 0.2.0 is being actively developed it is recommended to stay up-to-date with the latest code on the `master` branch.
  Passing `-e` to `pip3 install` will ensure that the `repESP` version installed reflects what's in your local git repo.
  To later update the code to latest version run: `git fetch --all ; git pull`.
* If you want to generate documentation, you'll need to append `[docs]` to the install command, i.e. run  `pip3 install .[docs]` instead.

# Contributing

`repESP` is free software released under the Gnu Public Licence version 3.
All contributions are very welcome, from improving the content of this overview to tweaking the code.

## If you cannot program ...

There's still plenty you can do!
First, have a look at GitHub's ["How to contribute"](https://guides.github.com/activities/contributing-to-open-source/#contributing).

### ... and you are familiar with `git`

* If you'd like to improve any of the documentation, fork this repository and submit a pull request.

### ... and you are not familiar with `git`

* If something doesn't work or should be improved, open an issue here on GitHub
* If you think a new feature would be nice, open an issue
* If you need a feature for your project (i.e. within a certain time frame), contact the maintainer directly to discuss whether it can be implemented time.
  Then we'll open an issue.

## If you can program

This program needs some work to get [where it wants to be](#where-is-repesp-going).
If you want to help out, please coordinate with the maintainer.
Developer instructions can be found [here](dev/README.md).
