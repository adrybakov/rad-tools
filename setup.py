# RAD-tools - Sandbox (mainly condense matter plotting).
# Copyright (C) 2022-2024  Andrey Rybakov
#
# e-mail: anry@uv.es, web: rad-tools.org
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from setuptools import find_packages, setup

from radtools import __version__

__version__ = __version__

with open("README.rst", "r") as file:
    long_description = ""
    for line in file:
        long_description += line

setup(
    name="rad-tools",
    version=__version__,
    description="Spin Hamiltonians, magnons and condense matter post-processing.",
    long_description=long_description,
    author="Andrey Rybakov",
    author_email="rybakov.ad@icloud.com",
    license="MIT license",
    url="https://rad-tools.org/en/stable/",
    download_url="https://github.com/adrybakov/rad-tools.git",
    packages=find_packages(),
    scripts=[
        "scripts/rad-plot-tb2j.py",
        "scripts/phonopy-plotter.py",
        "scripts/rad-extract-tb2j.py",
        "scripts/rad-make-template.py",
        "scripts/rad-plot-dos.py",
        "scripts/rad-plot-fatbands.py",
        "scripts/rad-identify-wannier-centres.py",
        "scripts/compute-energies.py",
        "scripts/rad-plot-tb2j-magnons.py",
    ],
    install_requires=["numpy", "matplotlib", "scipy", "tqdm", "termcolor"],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
    python_requires=">=3.8",
)
