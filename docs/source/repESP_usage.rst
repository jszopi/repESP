Usage notes
===========

The library is ready for use, but its API is not yet stable and there is room for improvement in some areas.
Known issues and considered enhancements are listed in the `Issues tab on GitHub`_; that is also the best place to ask questions about unexpected behaviour, unintuitive interfaces, unclear documentation, as well as put aspects of the library up for discussion.


External dependencies
---------------------

The library has a runtime dependency on the `resp` program and it operates on data produced by other programs.
See the `Prerequisites section`_ of the main GitHub page for information regarding obtaining and using these programs.

Runtime dependency
^^^^^^^^^^^^^^^^^^

Functions performing charge fitting, which are currently completely contained in the `resp_wrapper` module, invoke the `resp` program, so the program must be located in a directory listed in your `PATH` variable.
Functions in other modules operate on the output of the `resp` program.

Input data sources
^^^^^^^^^^^^^^^^^^

Functions which infer equivalence for use by the `resp` program require the output of the `antechamber` and `respgen` programs from the AmberTools suite of programs.
The Gaussian program is required to produce the ESP mesh and field, as well as the NPA charges.
These programs are not currently invoked by any of the library functions, so they do not need to be discoverable by your Python interpreter.


Units of measurement
--------------------

This library uses atomic units for all physical properties.
If input to your program is in other units, you should convert as early in your code as possible.
If you require output in other units, you should convert as late as possible.

If your existing code already uses different units or uses libraries that operate on units other than atomic, you must carefully alternate and convert between units.
It is easy for a mix-up to happen, resulting in difficult-to-spot numerical errors.
This library aims to prevent such a class of errors by requiring the units to be expressed explicitly in the code.
Instead of passing "naked" `float` values to this library, you should use the provided classes to represent physical properties.
For example, if you mean to specify a bond length of 0.5 a.u., use ``types.Dist(0.5)`` and if you mean to specify a charge of -0.1 e, use ``charges.Charge(-0.1)``.

The correct use of physical properties can then be enforced by the means of static type checking (see section `Type checking user code`_ below).
The type checker will complain if you try to pass a naked `float` or a `Dist` object into a library function expecting a charge::

    AtomWithCharge(6, 0.2)
    # error: Argument 2 to "AtomWithCharge" has incompatible type "float"; expected "Charge"
    AtomWithCharge(6, Dist(0.5))
    # error: Argument 2 to "AtomWithCharge" has incompatible type "Dist"; expected "Charge"

The types representing physical properties inherit from `float`, meaning that you will be able to use them in any context you can use a `float`.
Thus, the above `AtomWithCharge` construction executes just fine despite the type errors.
The types can also be passed into any other function accepting a `float`, for example::

    >>> numpy.mean([Dist(0.2), Dist(1.2)])
    0.7

Caveats
^^^^^^^

The type system for values representing physical properties can be considered an experiment and the implementation is a work in progress.
Some arithmetic operations may not result in the correct object type and will fall back onto the `float` type::

    Dist(0.5) + Charge(0.3)
    # No type checking error, addition is implemented for `float`s.

The fallback onto float may also result in false positives::

    c = Charge(0.5) - Charge(0.2)
    reveal_type(c)  # (this is a `mypy` command, it would be a runtime error)
    # error: Revealed type is 'builtins.float'
    AtomWithCharge(6, c)
    # error: Argument 2 to "AtomWithCharge" has incompatible type "float"; expected "Charge"

`Issue #49`_ tracks the shortcomings in the type system for physical properties.


User code correctness
---------------------

Precondition checking in Python
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Python is a language with a dynamic type system, where type checks tend to be limited to "duck typing".
The language allows to pass arguments of any type to a library function without providing any feedback at the point of making the function call.
One approach to ensuring correct use of libraries is for authors to write unambiguous and thorough documentation and for library users to abide by it perfectly.

Runtime precondition checks are another common approach, but implementing them at the API level is tedious.
Usually, precondition violations will manifests themselves deeper down the call stack, resulting in errors that are more difficult to interpret by library users.

Motivation for static typing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Statically-typed languages eliminate the need for a large class of precondition checks.
In Python, this approach can be realized with type annotations, which are typically checked before running user code.

The types expected and produced by this library can be found on every function call in the API documentation.
The library carries out no runtime checks of preconditions already expressed by the types; rather, the library may produce erroneous results if the input arguments do not follow the documented types.

Type checking user code
^^^^^^^^^^^^^^^^^^^^^^^

To avoid passing incorrect types, users are advised to check their code for type errors with a static type checker, like ``mypy``.
As an exercise, try introducing a type error into the example script and observe the output of::

    mypy --strict example.py

Ideally, you would use the ``--strict`` flag when checking your own code too.
However, if your code grows beyond a simple demo and you start creating functions, the ``--strict`` flag will report errors unless you also add type annotations to your code.
This requirement may be impractical for many library users, so in the future the library will offer checking the types at runtime instead, see `Issue #48`_.

Dataclasses
-----------

In any programming language, complex preconditions can be expressed with a custom data type.
This library makes heavy use of dataclasses, which are a Python construct to aggregate related pieces of data.
Generally, aggregate data types do not have any hidden state.
This is also the case for dataclasses used in this library with one crucial exception --- successful construction of an object asserts that the provided combination of input data is valid.
For example, the construction of an `Atom` will fail if the provided atomic number is 0 or 200.


.. _`Issues tab on GitHub`: https://github.com/jszopi/repESP/issues
.. _`Prerequisites section`: https://github.com/jszopi/repESP#prerequisites
.. _`Issue #48`: https://github.com/jszopi/repESP/issues/48
.. _`Issue #49`: https://github.com/jszopi/repESP/issues/49
