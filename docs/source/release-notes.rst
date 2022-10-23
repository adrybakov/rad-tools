*************
Release notes
*************

0.1
===
The big renaming passed.

0.1.0
-----
Scripts were renamed:

tb2j_plotter.py to :ref:`tb2j-plotter.py <tb2j-plotter>`

tb2j_prefractor.py to :ref:`tb2j-prefractor.py <tb2j-refractor>`

phonopy_plotter.py to :ref:`phonopy-plotter.py <phonopy-plotter>`

Modules were renamed:

file_logic to :doc:`model <api/rad_tools.exchange.model>`

template_logic to :doc:`template <api/rad_tools.exchange.template>`

map_logic to :doc:`map <api/rad_tools.map>`

tb2j_tools was renamed to :doc:`exchange <api/rad_tools.exchange>`

Module :doc:`map <api/rad_tools.map>` was moved out of 
:doc:`exchange <api/rad_tools.exchange>`.


0.0
===
Preliminary stage of the project, the main problem gere is a messy organisation.

0.0.3
-----
Add possibility to make draft of the template file form TB2J file in
:ref:`rad-make-template.py <rad-make-template>` script.

0.0.2
-----
Add :ref:`rad-make-template.py <rad-make-template>` script. 
Fix bugs in :ref:`tb2j-plotter.py <tb2j-plotter>`.

0.0.1
-----
Change versioning style, correct bugs in template logic.


0.0.0.10
--------
Add :ref:`tb2j-refractor.py <tb2j-refractor>` script.

0.0.0.9
-------
Better help messages in :ref:`tb2j-plotter.py <tb2j-plotter>` script.

0.0.0.8
-------
Add possibility to plot parameters vs distance from the center of the molecule
to the center of the bond (see 
:ref:`mode <tb2j-plotter_mode>` and 
:ref:`substrate-atoms <tb2j-plotter_substrate-atoms>`).

Add argument to :ref:`tb2j-plotter.py <tb2j-plotter>` for title for the pictures 
(see :ref:`title <tb2j-plotter_title>`).

0.0.0.7
-------
Add the :ref:`phonopy-plotter.py <phonopy-plotter>` script.

0.0.0.6
-------
Add arguments :ref:`scale-data <tb2j-plotter_scale-data>` and 
:ref:`scale-atoms <tb2j-plotter_scale-atoms>` to the 
:ref:`tb2j-plotter.py <tb2j-plotter>`.

0.0.0.5
-------
Fix the problem with the :py:mod:`exchange` docs. 

0.0.0.4
-------
First release with fully working documentation.
