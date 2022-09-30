from importlib.metadata import entry_points
from setuptools import setup, find_packages

__version__ = "0.0.0.1"

long_description = """Collection of scripts from my PhD."""

setup(
    name='rad-tools',
    version=__version__,
    description='Collection of scripts from my PhD',
    long_description=long_description,
    author='Andrey Rybakov',
    author_email='rybakov.ad@icloud.com',
    license='MIT license',
    url='https://github.com/adrybakov/rad-tools',
    download_url='https://github.com/adrybakov/rad-tools.git',
    packages=find_packages(),
    scripts=[
        'scripts/tb2j_plotter.py'
    ],
    install_requires=[
        'numpy', 'matplotlib'
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
    ],
    python_requires='>=3.6',
)
