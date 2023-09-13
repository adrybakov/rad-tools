******************
What is RAD-tools?
******************

The package covers post-processing scenarios for the 
results of |QE|_, |TB2J|_ and |Wannier90|_, as well as provides some custom scripts.

The main focus of the package is leaning towards the treatment of spin Hamiltonian 
with the parameters obtained form any source. It provides the means of solveng the 
magnon problem, simplify the notation changes, and aimed at the graphical as well as text output.

Besides the spin Hamiltonian it realises the general :ref:`Crystal <guide_crystal_crystal>`, 
:ref:`Lattice <guide_crystal_lattice>`, :ref:`Atom <guide_crystal_atom>` and :ref:`Kpoints <guide_crystal_kpoints>`
routines, which are essential for the treatment of the condense matter systems.
We provide automatic and unique choice of the cell for :ref:`Bravais lattices<library_bravais-lattices>`, 
which empowers the reproducibility of the results.

It is expected to be used in two ways:

* :ref:`As Python module <module-guide>`
* :ref:`Via scripts <scripts-guide>` (i.e. usage from command line)

.. note::
  The script`s output is colored, however we respect `NO_COLOR <https://no-color.org/>`_.
