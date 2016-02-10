# Implement research ideas

* Read in ESP fitting datapoints (.esp) produced by Gaussian for restricted ESP (rESP) fitting.
  Make them compatible with the original grid and check if values at each point match.
  Enable the same analysis for this subset of points.
  Calculate fit quality at such datapoints as well as the full grid to compare with Gaussian's printed RMS and RRMS in .log files.

## Analysis

# Code

## Features

### Distance from QTAIM atoms

Yield distance from QTAIM atom (QTAIM extraction only yields a field of atomic labels).
This might require splitting `_dist_func` but its design was intentional to prevent double execution and a lot of effort went into designing the program around it.
So it should preferably stay, even if its functionality will be partly redundant to calculating the atomic label field first and then the distance with a separate function.

## Testing

* Testing `cube_helpers` is not very modular.
    It can be broken down a bit, but maybe that should follow after breaking down the module itself, see the [refactoring](#refactoring) section.
    Separating the tests from their initialization is tricky, a custom test suite probably isn't the solution.
    Maybe the objects necessary for testing should be pickled?
    A separate test module could test if building the cube results in the same pickled object.

## Refactoring

* The graphs.plot function needs to be refactored to decouple filtering and allow focusing on the appearance of the graph more.

* Break down `cube_helpers` module.

* I'm slightly concerned about the objects being modified as the program runs.
    The concern arose when I was thinking about testing, so it may not be justified for normal runs of the program.

### Logical relations between objects / hierarchy

This will probably break tests and generally take long, so probably shouldn't be done before I have a some of my analysis ideas working.

A `Molecule` is one of the fundamental objects and it's a list of atoms.
`Atom`s have their charges, but the parent `Molecule`s should also contain information about which charges have been extracted.

A `Grid` is only extracted from a cube.

A `Field` is a combination of a `Molecule` and a `Grid`.
It can be extracted from a cube but is also a result of various calculations.
For a given `Molecule`, it should be possible to list all the derived `Field`s.

### Avoiding recalculations

This can be done either by keeping the results as attributes or memoizing function calls.
The best solution will vary on a case-by-case basis but should follow clearly from the established beforehand logical structure.
Maybe memoization should not check just the ID but some of the properties of object arguments?
