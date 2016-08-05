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
* 'compromise' charges

# How to use `repESP`

`repESP` implemented a lot of functionalities, which researchers can use to re-assess the various methods of assigning partial charges.
Some of these functionalities are exposed to users as command line (CLI) scripts, which are relatively easy to use.
The scripts perform well-defined, modular tasks and can be woven together to create one's own work flow (UNIX style).

Other functionalities are less user friendly and require modifying a few input values in Python scripts
The plan is to transform them into command line scripts as well.

The following chart presents the data flow that `repESP` can support.
The monospaced text above arrows shows the names of the relevant command line scripts.

![Data flow](https://github.com/jszopi/repESP/blob/master/data_flow.png)

Remember to run the scripts with the `.py` extension.
Start by running with the option `-h` to see a help message and usage instructions.

# Where is `repESP` going?

Every research project is different and the CLI scripts cannot cater for all but the most standard purposes.
Therefore, `repESP` will also attempt to provide an API/library, so that researchers can write their own Python scripts.
They will be able to take advantage of the general capabilities of `repESP`, e.g. reading in charges from Gaussian output and reproducing the ESP.
The high level of abstraction and clear documentation should make it usable by anyone comfortable with Python.
Ideally, the user will only have to code elements specific to their research project.

# Installation

## Prerequisites

`repESP` needs input from computational chemistry software to obtain the ESP and partial charges.
For the foreseeable future the only supported program is Gaussian 09 (adapting the program for older versions may be easy, please open an issue if you want this).

To run, `repESP` requires Python with certain packages and two other free programs.

### Python

* Python 3.5
* Packages: `scipy`, `matplotlib`, `fortranformatter` (install using `pip3.5`, which comes with your Python)

### `resp` 

`resp` is a free program necessary to perform ESP-fitting.
Download link can be found at the bottom of this [documentation page](http://upjv.q4md-forcefieldtools.org/RED/resp/).
To install you'll have to compile the Fortran source and put the resulting binary in your `PATH`.

### AmberTools

`resp` and `repESP` need information about atom equivalence in the form of `.respin` files.
AmberTools is a suite of free programs ([link](http://ambermd.org/AmberTools16-get.html)).
However, it requires compilation from C, C++ and Fortran sources, which proved tricky on MacOS.

To create the necessary input for a methane molecule, first create the `.ac` file from Gausian `.log` using `antechamber`:

`antechamber -i methane.log -fi gout -o methane.ac -fo ac`

Then run `respgen` twice:

`respgen -i methane.ac -o methane.respin1 -f resp1`
`respgen -i methane.ac -o methane.respin2 -f resp2`

## Install

The following instructions are for UNIX-like systems, like MacOS and Linux.

* Choose and go to installation directory (e.g. `~/bin`)
* Get the code: `git clone https://github.com/jszopi/repESP`
* From inside the repESP directory run: `python3.5 setup.py develop`
* Make the scripts callable: `export PATH=$PATH:your_scripts_path`

To update the code to latest version run: `git fetch --all ; git pull`

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
The most burning need is defining the API and filling the discovered gaps.
This will allow to release a stable version of `repESP`.

The remaining kludgey scripts need to be rewritten as standalone CLI scripts.
This is crucial to make the program less user-hostile.
