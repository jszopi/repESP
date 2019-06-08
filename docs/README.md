# API documentation

The files in the `source` directory were originally generated with:

```sh
sphinx-apidoc -fMe -o source ../repESP
```

However, it's not recommended to run the command again, as it will overwrite the modifications made to these files.

# Adding modules

When a new file is added to the library, a new corresponding template file must be added.
Use one of the existing template files, e.g. `source/repESP.types.rst`, which contents are currently:

```
repESP.types module
===================

.. automodule:: repESP.types
    :members:
    :undoc-members:
    :show-inheritance:
```

You will also need to add the module to the `toctree` in `source/repESP_modules.rst`.
