In scripts use imports from the package (rather than relative), e.g.:

	from repESP import cube_helpers
	from repESP.graphs import plot

Such scripts can be run from any directory but the package needs to be installed (`setup.py`) beforehand.
For development, run only once

	sudo python3 setup.py develop

to link the files in the `repESP` directory and hence avoid reinstalling the package after every change.
[Source](http://tjelvarolsson.com/blog/begginers-guide-creating-clean-python-development-environments/).
