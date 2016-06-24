# Code based on LPTHW:
# http://learnpythonthehardway.org/book/ex46.html

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'name': 'repESP',
    'version': '0.1.0',
    'description': 'ESP field reproduction from partial charges and analysis '
                   'of resulting discrepancies',
    'author': 'Jan Szopinski',
    'author_email': 'jszopi@users.noreply.github.com',
    'url': 'https://github.com/jszopi/repESP',
    'packages': ['repESP'],
    'install_requires': ['scipy', 'matplotlib', 'fortranformat'],
    'license': 'GPLv3',
}

# TODO Additional files will need to be added when demo scripts are created
# https://docs.python.org/3/distutils/setupscript.html#installing-additional-files

setup(**config)
