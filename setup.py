"""
Module computing Resource Balance Analysis.
This file is adapted from

https://github.com/pypa/sampleproject/
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

with open(path.join(here, 'rba', '_version.py'), 'r') as f:
    version = f.readline().split("'")[1]

with open(path.join(here, 'rba', '_authors.py'), 'r') as f:
    authors = f.readline().split("'")[1]

setup(
    name='rbapy',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=version,

    description='Package for automated generation of bacterial Resource Balance Analysis (RBA) models and simulation of RBA models',
    long_description=long_description,

    # The project's main homepage.
    url='https://sysbioinra.github.io/RBApy/',

    # Author details
    author=authors,
    author_email='anne.goelzer@inra.fr',

    # Choose your license
    license='GPLv3',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
    ],

    # What does your project relate to?
    keywords=[
        'metabolism',
        'Resource Balance Analysis',
        'molecular biology',
        'cell biology',
        'biochemistry',
        'systems biology',
        'computational biology',
        'mathematical modeling',
        'numerical simulation',
    ],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),

    # Alternatively, if you want to distribute just a my_module.py, uncomment
    # this:
    #   py_modules=["my_module"],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=[
        'biopython',
        'lxml',
        'numpy',
        'pandas',
        'python-libsbml',
        'scipy',
        'requests'
    ],

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    extras_require={
        'cplex': ['cplex'],
        'cplex_conv_opt': ['conv_opt[cplex]'],
        'cplex_optlang': ['optlang', 'cplex'],

        'glpk': ['glpk'],
        'swiglpk': ['swiglpk'],
        'glpk_conv_opt': ['conv_opt'],
        'glpk_optlang': ['optlang'],

        'gurobi': ['gurobipy'],
        'gurobi_conv_opt': ['conv_opt[gurobi]'],
        'gurobi_optlang': ['optlang', 'gurobipy'],

        'dev': ['check-manifest'],
        'test': ['pytest', 'pytest-cov', 'coverage'],
    },

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={
        #'sample': ['package_data.dat'],
    },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    #data_files=[('my_data', ['data/data_file'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'generate-rba-model=rba.cli.generate_rba_model:main',
            'solve-rba-model=rba.cli.solve_rba_model:main',
        ],
    },
)
