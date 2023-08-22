*******************************************
Tools for the package setup and development
*******************************************

All scripts, which are intended to be use in the development/deployment process
are located here, unless they belong in the docs/ or tests/ directories.

Makefile is not located here, but heavily relies on the scripts in this directory.

All scripts are intended to be run from the root directory of the project.


Tools
=====

generate-scripts-docs.py
------------------------

Script generate Arguments section in the documentation and create_parser() function
in the script implementation based on the manager() function in the script implementation.
this script.

new-script.py
-------------

Create template files for the new script.

plot-bravais-lattice.py
-----------------------

Plot all examples of the bravais lattice.

plot-data-structure.py
-----------------------

Plot the data structure of the package. 

plot-notation.py
----------------

Plot notation tree picture.

plot-script-guide.py
--------------------

Plot pictures from the script user-guide (examples folder).

prepare-release.py
------------------

Prepare the release of the package. Never release without the successful run of