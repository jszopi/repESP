# Code based on LPTHW:
# http://learnpythonthehardway.org/book/ex46.html

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

docs_require = ['sphinx']
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
        'console_scripts': [
            'adjusted_charges=scripts.adjusted_charges:main',
            'average=scripts.average:main',
            'cavity=scripts.cavity:main',
            'dipole=scripts.dipole:main',
            'eval_fit=scripts.eval_fit:main',
            'field_diff=scripts.field_diff:main',
            'fit_dependence=scripts.fit_dependence:main',
            'fit_points=scripts.fit_points:main',
            'flexibility=scripts.flexibility:main',
            'plot_fit_dependence1=scripts.plot_fit_dependence1:main',
            'plot_fit_dependence2=scripts.plot_fit_dependence2:main',
            'rep_esp=scripts.rep_esp:main',
            'run_two-stage_resp=scripts.run_two_stage_resp:main',
        ]
    }
}

# TODO Additional files will need to be added when demo scripts are created
# https://docs.python.org/3/distutils/setupscript.html#installing-additional-files

setup(**config)
