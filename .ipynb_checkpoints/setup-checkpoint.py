#import setuptools
from __future__ import print_function
from numpy.distutils.core import setup, Extension

setup(name = "H2model",
      version = "0.1",
      description = "Toy model calculator",
      author = "Antonio Siciliano Lorenzo Monacelli",
      packages = ["H2model"],
      package_dir = {"H2model": "Modules"},
      license = "GPLv3")


def readme():
    with open("README.md") as f:
        return f.read()