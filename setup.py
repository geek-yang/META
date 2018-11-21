#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

with open('README.rst') as readme_file:
    readme = readme_file.read()

setup(
    name='meta',
    version='0.1.0',
    description="Meridional Energy Transport Analyzer, in short as META, is a python library for the calculation of meridional energy transport and further analysis based on climate data. It is able to deal with the calculations of MET in both the atmosphere (AMET) and ocean (OMET). In addition, it provides diagnostic modules to perform statistical operation on MET.",
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
