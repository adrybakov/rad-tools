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
    description="Scripts for condense matter calculations and post-processing.",
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
        "License :: OSI Approved :: MIT License",
    ],
    python_requires=">=3.8",
)
