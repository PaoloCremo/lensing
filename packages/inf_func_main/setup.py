# example:
# https://github.com/PyFstat/PyFstat/blob/master/setup.py

import setuptools

setuptools.setup(
    name="inf_func",
    packages=setuptools.find_packages(),
    install_requires=[
        #'',
        'numpy',
        # 'bilby',
        # 'bilby_pipe',
        'scipy',
        'mpmath',
        # '...',
        ]
)

import re
VERSIONFILE="inf_func/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))
~                                                                               