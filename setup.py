#!/usr/bin/env python3
from sys import version_info, stderr
from setuptools import setup, find_packages

NAME = 'pancat'
CURRENT_PYTHON = version_info[:2]
REQUIRED_PYTHON = (3, 10)

if CURRENT_PYTHON < REQUIRED_PYTHON:
    stderr.write(
        f"{NAME} requires Python 3.10 or higher and your current version is {CURRENT_PYTHON}.")
    exit(1)

setup(
    name=NAME,
    version='0.3.3',
    description='Tools for manipulating and visualizing GFA file format',
    url='https://github.com/dubssieg/pancat',
    author='dubssieg',
    author_email='dubois.siegfried@inria.fr',
    packages=find_packages(),
    package_data={'': ['template.html']},
    include_package_data=True,
    zip_safe=False,
    license="LICENSE",
    long_description=open("README.md", encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    install_requires=[
        'matplotlib',
        'seaborn',
        'numpy',
        'pyvis',
        'networkx',
        'mycolorpy',
        'levenshtein',
        'gfagraphs',
        'BubbleGun',
        'statsmodels',
        'resource',
        'Bio',
        'rich'
    ],
    entry_points={'console_scripts': ['pancat=pancat.main:main']}
)
