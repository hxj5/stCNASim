# setup.py


import sys
from codecs import open
from os import path
from setuptools import setup, find_packages



here = path.abspath(path.dirname(__file__))


# load configures
exec(open("./sccnasim/app.py").read())


# Get the long description from the relevant file
with open(path.join(here, "README.rst"), encoding='utf-8') as f:
    long_description = f.read()

    
# check Python version.
# currently, we tested Python 3.7 (not compatible), 3.11 (compatible).
# Ref: https://stackoverflow.com/questions/13924931/setup-py-restrict-the-allowable-version-of-the-python-interpreter
if sys.version_info[0] != 3 or sys.version_info[1] < 11:
    sys.exit("Sorry, only Python 3.11 or above is supported.")

    
# Alternatively, use `poetry` with `pyproject.toml` file to specify the 
# Python version, as `XClone` did.


# set "numpy<2"?
reqs = ["anndata", "intervaltree", "matplotlib", "numpy", "pandas", "pysam",
        "scipy", "seaborn", "statsmodels"]


setup(
    name = "sccnasim",

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version = VERSION,

    description = "scCNASim - Haplotype-aware simulation of somatic CNAs from single-cell and spatial transcriptomics",
    long_description = long_description,
    long_description_content_type = "text/markdown",

    # The project's main homepage.
    url = "https://github.com/hxj5/scCNASim",

    # Author details
    author = 'Xianjie Huang',
    author_email = 'xianjie5@connect.hku.hk',

    # Choose your license
    license='Apache-2.0',

    # What does your project relate to?
    keywords = ['CNA', 'Simulator', "single-cell", "spatial transcriptomics"],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages = find_packages(),

    # or use
    #packages = ["sccnasim"],
    #package_dir = {"": "."},

    #entry_points={
    #    'console_scripts': [
    #        'sccnasim = sccnasim.main:main'
    #    ],
    #},

    # Check Python version. Only pip >= 9.0.1 supports.
    # Ref: https://stackoverflow.com/questions/42238484/prevent-package-from-being-installed-on-old-python-versions/42792413
    python_requires = ">=3.11",

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    
    install_requires = reqs,

    py_modules = ['sccnasim']

    # buid the distribution: python setup.py sdist
    # upload to pypi: twine upload dist/...
)
