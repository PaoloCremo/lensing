# example:
# https://github.com/PyFstat/PyFstat/blob/master/setup.py

exec(open('inf_func/_version.py').read())

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
