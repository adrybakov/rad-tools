*************
Release notes
*************

0.2
===
Add script :ref:`rad-dos-plotter.py <rad-dos-plotter>`

0.2.12
------
Bugfix in :ref:`tb2j-plotter.py <tb2j-plotter>` script. 

0.2.11
------
Rename --template option in 
:ref:`tb2j-extractor.py <tb2j-extractor>` script into 
:ref:`--template-file <tb2j-extractor_template-file>`. 

0.2.10
------
Rename tb2j-refractor.py to :ref:`tb2j-extractor.py <tb2j-extractor>`

0.2.9
-----
Add pythonic :ref:`--verbose <tb2j-extractor_verbose-ref>` output
to the :ref:`tb2j-refractor.py <tb2j-extractor>` script.

0.2.8
-----
Add :ref:`--verbose <tb2j-extractor_verbose>` flag
to the :ref:`tb2j-refractor.py <tb2j-extractor>` script.

0.2.7
-----
Change the behaviour of the drawing in 
:ref:`tb2j-plotter.py <tb2j-plotter>` script.

0.2.6
-----
Add :ref:`-dmi <tb2j-extractor_dmi>` flag
to :ref:`tb2j-refractor.py <tb2j-extractor>` script.

0.2.5
-----
Add \|DMI\| as an output data type to :ref:`tb2j-plotter.py <tb2j-plotter>` 
(see: :ref:`--what-to-plot <tb2j-plotter_what-to-plot>`)

0.2.4
-----
Bug fix in :ref:`tb2j-refractor.py <tb2j-extractor>`. 
Only last exchange from template was shown, now all of them are.

0.2.3
-----
Change output behaviour in :ref:`tb2j-refractor.py <tb2j-extractor>`.
Now by default output is passed to the standart output stream.


0.2.2
-----
Add interactive mode to the :ref:`rad-dos-plotter.py <rad-dos-plotter>`.

0.2.1
-----

Correct output file name in :ref:`rad-dos-plotter.py <rad-dos-plotter>`.

0.1
===
The big renaming passed.

0.1.1
-----
Fix bugs in :ref:`tb2j-refractor.py <tb2j-extractor>`.

0.1.0
-----
Scripts were renamed:

tb2j_plotter.py to :ref:`tb2j-plotter.py <tb2j-plotter>`

tb2j_refractor.py to :ref:`tb2j-refractor.py <tb2j-extractor>`

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
Preliminary stage of the project, the main problem here is a messy organisation.

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
Add :ref:`tb2j-refractor.py <tb2j-extractor>` script.

0.0.0.9
-------
Better help messages in :ref:`tb2j-plotter.py <tb2j-plotter>` script.

0.0.0.8
-------
Add possibility to plot parameters vs distance from the center of the molecule
to the center of the bond (see 
:ref:`--mode <tb2j-plotter_mode>` and 
:ref:`--atoms <tb2j-plotter_atoms>`).

Add argument to :ref:`tb2j-plotter.py <tb2j-plotter>` for title for the pictures 
(see :ref:`--title <tb2j-plotter_title>`).

0.0.0.7
-------
Add the :ref:`phonopy-plotter.py <phonopy-plotter>` script.

0.0.0.6
-------
Add arguments :ref:`--scale-data <tb2j-plotter_scale-data>` and 
:ref:`--scale-atoms <tb2j-plotter_scale-atoms>` to the 
:ref:`tb2j-plotter.py <tb2j-plotter>`.

0.0.0.5
-------
Fix the problem with the :py:mod:`.exchange` docs. 

0.0.0.4
-------
First release with fully working documentation.
