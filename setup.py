from setuptools import setup, find_packages
from rad_tools import __version__

__version__ = __version__

with open('README.rst', 'r') as file:
    long_description = ''
    for line in file:
        long_description += line

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
        'scripts/tb2j-plotter.py',
        'scripts/phonopy-plotter.py',
        'scripts/tb2j-refractor.py',
        'scripts/rad-make-template.py',
        'scripts/rad-jarvis.py',
        'scripts/rad-dos-plotter.py'
    ],
    install_requires=[
        'numpy', 'matplotlib'
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
    ],
    python_requires='>=3.6',
)
