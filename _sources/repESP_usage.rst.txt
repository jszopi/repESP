Usage notes
===========

Caveats
-------

* This library uses atomic units for all physical properties. If you're dealing with
  input in other units, you should convert as early in your program as possible.
  If you require output in other units, you should convert as late as possible.
* This library attempts to enforce consistent use of values describing physical
  properties. You shouldn't need to use "naked" `float` values to represent
  physical properties and should use the provided classes instead, which are
  lightweight and inherit from `float`. For example, if you mean to specify a
  bond length of 0.5 a.u., you should use ``types.Dist(0.5)`` and if you mean
  to specify a charge of -0.1 e, you should use ``charges.Charge(-0.1)``. This
  aims to ensure that you never try to subtract distance from a charge or try
  to add values in different units. This area of the library may not be
  bulletproof and some arithmetic operations may not result in the correct
  classes and will fall back onto the `float` type, possibly causing type
  checker errors. Note that you will be able to use the inheriting types like
  `Dist` in any context you can use a `float`.
* Programs using this library are required to use correct types for input
  arguments. The types are documented here and additionally the user should run
  a static type checker, like ``mypy`` on their script beforehand. To practice,
  try introducing a type error into the example script and observe the output of::

      mypy --strict example.py

  Using a static type checker will save the user debugging tricky runtime
  errors which tend to occur later in the program. Even if the program execution
  is not stopped by a runtime error, the library may produce erroneous results
  if the passed-in arguments don't follow the documented types.
