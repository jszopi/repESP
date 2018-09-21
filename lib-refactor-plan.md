# Priorities

* Thought-through API.
  I haven't been spending much time on thinking through the recent decisions (NumericField being an example) even though these are admittedly of lower impact than the better-thought-through groundwork.
  The beauty or lack thereof of the API will only be clear once docs are complete and it's mini-battle tested in the scripts.
* Documentation
* No loss of functionality - it will be necessary to maintain two version in parallel, with some functionality only available through scripts.

## Primary

Required for the basic functionality in API.

* qout format, needed for resp but also fair format for charges
* resp
* flexibility

## Secondary

These are not crucial for the early version of the API. They are required to make the two scripts work, `field_diff` and `fit_points` (and possibly a bit of graphs needed for `plot_fit_dependence`). These are important scripts, so unless you can have to versions of the program running (or subsum 0.1 into 0.2 privately only for the benefit of those scripts, which is not elegant but will work).

* field relative difference - only need to implement field division
* field filtering
* graphs - may actually be decoupled strongly enough for this to be easy.

## Bits

These are quick but of lower priority. Some are required for other scripts: `cavity` and `dipole`.

* AIM parsers
* distance transform (`cavity`)
* dipole (and quadrupole) calculation (`dipole`)
* checking for NANs - how can they get introduced except for my calculations?
* Printing in various units (especially distances in Angstrom)

## External uses

* Add dev dependencies (mypy, snapshot testing, test coverage)
* Unless version 0.2.0 is ready for October (including all scripts), users should be able to run both version in parallel, probably requiring a virtualenv.
* If virtualenv is required (or I have the time for improvements), pipenv is probably better.

## Automation

* AC tests for scripts? Scripts may be easy enough to be done quickly at the end but should probably do one earlier to find out for sure.
* Docs - not compulsory, will prob start documenting later but should comment/stub now
* Dev checks are quite good, setting up GitHub CI would allow to reap the benefits.
