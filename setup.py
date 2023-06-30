#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

# access the version wihtout importing the SAPiN package
with open('SAPiN/version.py') as f: exec(f.read())

setup(
    name='SAPiN',
    version=__version__,

    description='summarize BAM read alignment by pileup or reads at each position in a tabulated way',

    url='https://github.com/Juke34/SAPiN',
    download_url='https://github.com/Juke34/SAPiN/archive/v' + __version__ +'.tar.gz',
    author='Jacques Dainat',

    license='GPL-3.0',
    packages=find_packages(),

    install_requires=[
        'gffutils >= 0.11.1', 'numpy>=1.22', 'matplotlib',
        ],
    include_package_data=True,

    entry_points={
        'console_scripts': ['sapin = SAPiN:main',
        ],
    }
)
