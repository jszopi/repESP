# Code based on LPTHW:
# http://learnpythonthehardway.org/book/ex46.html

from scripts.docs.script_list import scripts_to_hooks

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


def get_console_scripts(script_to_hook_dict):
    return [f"{script}={hook}" for script, hook in script_to_hook_dict.items()]


docs_require = ['sphinx', 'sphinx_rtd_theme']
test_requires = ['mypy', 'nose']
dev_requires = docs_require + test_requires

config = {
    'name': 'repESP',
    'version': '0.2.0-dev',
    'description': 'Reproducing molecular ESP from partial charges and more',
    'author': 'Jan Szopinski',
    'author_email': 'jszopi@users.noreply.github.com',
    'url': 'https://github.com/jszopi/repESP',
    'packages': ['repESP'],
    'install_requires': ['scipy', 'matplotlib', 'fortranformat', 'pandas'],
    'extras_require': {
        'docs': docs_require,
        'dev': dev_requires,
        'test': test_requires,
    },
    'python_requires': '>=3.7.0',
    'license': 'GPLv3',
    'entry_points': {
        'console_scripts': get_console_scripts(scripts_to_hooks),
    }
}

# TODO Additional files will need to be added when demo scripts are created
# https://docs.python.org/3/distutils/setupscript.html#installing-additional-files

setup(**config)
