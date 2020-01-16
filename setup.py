#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

with open('README.rst') as readme_file:
    readme = readme_file.read()

setup(
    name='meta',
    version='0.1.2',
    description="Meridional Energy Transport Analyzer, in short as META, is a python library designed for the calculation of meridional energy transport (MET) and extended analysis with climatological data sets. It is able to deal with the calculations of MET in both the atmosphere (AMET) and ocean (OMET). In addition, it provides diagnostic modules to perform statistical operations on MET.",
    long_description=readme + '\n\n',
    author="Yang Liu",
    author_email='y.liu@esciencecenter.nl',
    url='https://github.com/geek-yang/meta',
    packages=[
        'meta',
    ],
    package_dir={'meta':
                 'meta'},
    include_package_data=True,
    install_requires=[
        'ncl>=6.4.0'
    ]
    license="Apache Software License 2.0",
    zip_safe=False,
    keywords='meta',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    test_suite='tests',
    setup_requires=[
        # dependency for `python setup.py test`
        'pytest-runner',
        # dependencies for `python setup.py build_sphinx`
        'sphinx',
        'recommonmark'
    ],
    tests_require=[
        'pytest',
        'pytest-cov',
        'pycodestyle',
    ],
    extras_require={
        'dev':  ['prospector[with_pyroma]', 'yapf', 'isort'],
    }
)
