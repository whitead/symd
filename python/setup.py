import os
from glob import glob
from setuptools import setup

exec(open("symd/version.py").read())

with open("../README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="symd",
    version=__version__,
    description="Symmetric Molecular Dynamics",
    author="Andrew White",
    author_email="andrew.white@rochester.edu",
    url="https://whitead.github.io/symd/",
    license="MIT",
    packages=["symd", "symd.data"],
    install_requires=[
        "selfies >= 2.0.0",
        "numpy",
        "requests",
        "tqdm",
        "ratelimit",
        "rdkit-pypi",
        "scikit-learn",
        "skunk >= 0.4.0",
    ],
    test_suite="tests",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
)
