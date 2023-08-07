.. _guide_exchange:

********
Exchange
********

For the full reference see :ref:`api_exchange`.

.. currentmodule:: radtools

Exchange module describes the spin Hamiltonian and is build around of the 
:py:class:`.SpinHamiltonian` class. It stores bonds between :py:class:`.Atom` 
:math:`i`, which is located in :math:`(0, 0, 0)` unit cell and :py:class:`.Atom` :math:`j`,
which is located in :math:`(i, j, k)` unit cell and corresponding exchange parameters.

The notation of the exchange Hamiltonian is an important issue, since without clear notation 
exchange parameters does not make much sense. :py:class:`.SpinHamiltonian` can support any 
notation out of the most common ones. Detailed description of the notation is given in 
correspondent section: :ref:`guide_exchange_notation`. 
We encourage you to read it once with full attention.

The main building block of the exchange Hamiltonian is exchange parameter 
(:math:`\boldsymbol{J}`). It is implemented in a separate class :py:class:`.ExchangeParameter`.
For the guide on exchange parameter see:

.. toctree::
    :maxdepth: 1
    
    parameter


.. _guide_exchange_notation:

Notation
========

The notation, which is considered to be the standard for the RAD-tools:

.. math::

    H = -\sum_{i,j} \hat{\boldsymbol{S}}_i \cdot \boldsymbol{J}_{i,j} \hat{\boldsymbol{S}}_j

It has double-counting, spins are not normalized, and the exchange parameter is
positive for ferromagnetic order.

.. note::
    The "standard" here does not mean that the :py:class:`.SpinHamiltonian` is
    always has this notation. For example when :py:class:`.SpinHamiltonian` is
    read from |TB2J|_ file (:py:func:`.read_tb2j_model`) it has the notation of TB2J.

There are five properties, that has to be defined in order to describe the 
exchange Hamiltonian`s notation:

* :py:attr:`.double_counting`
    Whether both pairs :math:`(i, j)` and :math:`(j, i)` are included in the sum. 

    Double counting is avoided it is usually indicated under the sum sign:

    .. math::

        \sum_{i < j} \text{ or } \sum_{i > j}

    When it is not present there are no indication in the sum: 

    .. math::

        \sum_{ij}

    However, some authors are using 
    
    .. math::
        
        \sum_{<ij>} 
    
    as an indication of the avoided double counting, which we find confusing 
    and discourage you to use it in this way. 
    In textbooks this indication  is usually imply that only near-neighbors
    are included in the sum, which is not the same as avoiding double counting 
    (consecutively :math:`\sum_{<<ij>>}` is used to indicate next-nearest neighbors and so on).

    .. note::

        Indication

        .. math::

            \sum_{i \ne j}

        Does not mean that double counting is avoided. 
        It specifies that there is no exchange between the same atom.

    
* :py:attr:`.spin_normalized`
    :math:`\boldsymbol{S}_i` in the exchange Hamiltonian is viewed as spin operator 
    or as classical spin vector. In the latter case it could be normalized to 1 or not. 
    In case of normalisation exchange parameter absorbs the factor 
    :math:`\vert\boldsymbol{S}_i\vert \vert\boldsymbol{S}_j\vert`.
* :py:attr:`.factor_one_half`
    Factor :math:`1/2` is sometimes included in order to correct for double counting. If 
    it is included, then it is simply written before the sum sign. It is usually used
    when double counting is present.
* :py:attr:`.factor_two`
    Factor :math:`2` is sometimes included in order to account for double counting. If 
    it is included, then it is simply written before the sum sign. It is usually used
    when double counting is avoided.
* :py:attr:`.minus_sign`
    Whether the minus sign is included in the Hamiltonian. 

.. caution::

    We would like to note that not all authors are thoughtful with the definition
    of the Hamiltonian, so additional care is required when reading the literature.

RAD-tools utilize those five properties in order to define the notation of the
exchange Hamiltonian. During the creation of the SpinHamiltonian object the notation 
is deliberately not defined, because it depends on your interpretation. Therefore,
the notation has to be defined explicitly by you. If the notation is not defined, 
and the you are trying to use the properties and methods, which expect the notation
to be defined, then the :py:exc:`.NotationError` is raised.

The notation could be defined in two ways:

* By setting each individual property. For example:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> hamiltonian.double_counting = True
    >>> hamiltonian.spin_normalized = False
    >>> hamiltonian.factor_one_half = False
    >>> hamiltonian.factor_two = True
    >>> hamiltonian.minus_sign = True

* By setting the :py:attr:`.notation` property directly. For example:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> hamiltonian.notation = (True, False, False, True, True)
    >>> hamiltonian.notation
    H = -2 sum_{i,j} S_i J_ij S_j
    Double counting is present.
    Spin vectors are not normalized.
    (True, False, False, True, True)
    >>> hamiltonian.double_counting
    True
    >>> hamiltonian.spin_normalized
    False
    >>> hamiltonian.factor_one_half
    False
    >>> hamiltonian.factor_two
    True
    >>> hamiltonian.minus_sign
    True

Once the notation or any of the individual properties are set, 
the following redefinition of the notation or corresponding property will change exchange
parameters. See :ref:`Examples <guide_exchange_notation-examples>`.

If you want to change the notation once is set, but keep the parameters intact use 
:py:meth:`.SpinHamiltonian.set_interpretation` method:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> hamiltonian.double_counting = True
    >>> hamiltonian.spin_normalized = False
    >>> hamiltonian.factor_one_half = False
    >>> hamiltonian.factor_two = True
    >>> hamiltonian.minus_sign = True
    >>> hamiltonian.notation
    H = -2 sum_{i,j} S_i J_ij S_j
    Double counting is present.
    Spin vectors are not normalized.
    (True, False, False, True, True)
    >>> hamiltonian.set_interpretation(double_counting=False, spin_normalized=True, factor_one_half=False, factor_two=False, minus_sign=True)
    >>> hamiltonian.notation
    H = -sum_{i>=j} S_i J_ij S_j
    No double counting.
    Spin vectors are normalized to 1.
    (False, True, False, False, True)

Predefined notations
--------------------

Order: (double counting, spin normalized, factor 1/2, factor 2, minus sign).

* Standard
    (True, False, False, False, True)

    .. math::
        H = -\sum_{i,j} \hat{\boldsymbol{S}}_i \cdot \boldsymbol{J}_{i,j} \hat{\boldsymbol{S}}_j

    where double counting is present (:math:`ij` and :math:`ji` is in the sum).
    Spin vectors are **not** normalized.
* |TB2J|_
    (True, True, False, False, True)

    .. math::
        H = -\sum_{i,j} \hat{\boldsymbol{S}}_i \cdot \boldsymbol{J}_{i,j} \hat{\boldsymbol{S}}_j

    where double counting is present (:math:`ij` and :math:`ji` is in the sum).
    Spin vectors are normalized to 1.
* SpinW
    (True, False, False, False, False)

    .. math::
        H = \sum_{i,j} \hat{\boldsymbol{S}}_i \cdot \boldsymbol{J}_{i,j} \hat{\boldsymbol{S}}_j

    where double counting is present (:math:`ij` and :math:`ji` is in the sum).
    Spin vectors are **not** normalized.

They could be set directly through the :py:attr:`.notation` property:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> hamiltonian.notation = 'standard'
    >>> hamiltonian.notation
    H = -sum_{i,j} S_i J_ij S_j
    Double counting is present.
    Spin vectors are not normalized.
    (True, False, False, False, True)
    >>> hamiltonian.notation = 'tb2j'
    >>> hamiltonian.notation
    H = -sum_{i,j} S_i J_ij S_j
    Double counting is present.
    Spin vectors are normalized to 1.
    (True, True, False, False, True)
    >>> hamiltonian.notation = 'spinw'
    >>> hamiltonian.notation
    H = sum_{i,j} S_i J_ij S_j
    Double counting is present.
    Spin vectors are not normalized.
    (True, False, False, False, False)

.. _guide_exchange_notation-examples:

Examples
--------
Setting one of the predefined notations:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> hamiltonian.notation = "standard"
    >>> hamiltonian.notation
    H = -sum_{i,j} S_i J_ij S_j
    Double counting is present.
    Spin vectors are not normalized.
    (True, False, False, False, True)
    >>> hamiltonian.notation = "TB2J"
    >>> hamiltonian.notation
    H = -sum_{i,j} S_i J_ij S_j
    Double counting is present.
    Spin vectors are normalized to 1.
    (True, True, False, False, True)
    >>> hamiltonian.notation = "SpinW"
    >>> hamiltonian.notation
    H = sum_{i,j} S_i J_ij S_j
    Double counting is present.
    Spin vectors are not normalized.
    (True, False, False, False, False)

Setting the notation:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> Cr = rad.Atom("Cr", spin=1.5)
    >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=1), Cr, Cr, (1, 0, 0))
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    1.0
    >>> # For the first time interpretation is set,
    >>> # values of exchange are not changed
    >>> hamiltonian.notation = "standard"
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    1.0
    >>> # Once the notation is set the values
    >>> # are changing if the notation is changed again.
    >>> hamiltonian.notation = "TB2J"
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    2.25
    >>> hamiltonian.notation = "SpinW"
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    -1.0
    >>> hamiltonian.notation = "standard"
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    1.0

Setting individual properties:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> Cr = rad.Atom("Cr", spin=1.5)
    >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=1), Cr, Cr, (1, 0, 0))
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    1.0
    >>> # For the first time interpretation is set,
    >>> # values of exchange are not changed
    >>> hamiltonian.minus_sign = True
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    1.0
    >>> # Once the property is set the values
    >>> # are changing if the property is changed again.
    >>> hamiltonian.minus_sign = False
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    -1.0

Changing individual properties:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> Cr = rad.Atom("Cr", spin=1.5)
    >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=1), Cr, Cr, (1, 0, 0))
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    1.0
    >>> hamiltonian.notation = "standard"
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    1.0
    >>> hamiltonian.double_counting, hamiltonian.spin_normalized, hamiltonian.factor_one_half, hamiltonian.factor_two, hamiltonian.minus_sign
    (True, False, False, False, True)
    >>> hamiltonian.minus_sign = False
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    -1.0
    >>> hamiltonian.factor_one_half, hamiltonian.factor_two
    (False, False)
    >>> hamiltonian.factor_one_half = True
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    -2.0
    >>> hamiltonian.factor_one_half, hamiltonian.factor_two
    (True, False)
    >>> hamiltonian.factor_two = True
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    -1.0
    >>> # Note that the values are switched to False,
    >>> # since factor one half and two are cancelling each other
    >>> hamiltonian.factor_one_half, hamiltonian.factor_two
    (False, False)
    >>> hamiltonian.spin_normalized = True
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    -2.25
    >>> hamiltonian.double_counting = False
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    -4.5

Crystal of the exchange Hamiltonian
===================================

Exchange Hamiltonian need the crystal structure to be defined. 
:py:class:`.SpinHamiltonian` is a child of :py:class:`.Crystal` and inherits all its
properties and methods. Thus, for the creation of the crystal structure we refer you to the 
:ref:`guide_crystal_crystal` and :ref:`guide_crystal_atom` guides. 
Any property, which is related to the structure is expected to be called directly from the 
:py:class:`.SpinHamiltonian` instance. For example:

.. doctest::

    >>> import radtools as rad
    >>> lattice = rad.lattice_example("TET")
    >>> crystal = rad.Crystal(lattice)
    >>> hamiltonian = rad.SpinHamiltonian(crystal)
    >>> hamiltonian.variation
    'TET'
    >>> hamiltonian.a1
    array([3.14159265, 0.        , 0.        ])
    >>> hamiltonian.cell
    array([[3.14159265, 0.        , 0.        ],
           [0.        , 3.14159265, 0.        ],
           [0.        , 0.        , 4.71238898]])

The crystal of the Exchange Hamiltonian can be access through the 
:py:attr:`.SpinHamiltonian.crystal` attribute. 
Note, that it returns an independent instance of the :py:class:`.Crystal` class:

.. doctest::

    >>> import radtools as rad
    >>> lattice = rad.lattice_example("TET")
    >>> crystal = rad.Crystal(lattice)
    >>> hamiltonian = rad.SpinHamiltonian(crystal)
    >>> crystal = hamiltonian.crystal
    >>> crystal.cell
    array([[3.14159265, 0.        , 0.        ],
           [0.        , 3.14159265, 0.        ],
           [0.        , 0.        , 4.71238898]])
    >>> hamiltonian.cell
    array([[3.14159265, 0.        , 0.        ],
           [0.        , 3.14159265, 0.        ],
           [0.        , 0.        , 4.71238898]])
    >>> crystal.cell = [[1,0,0],[0,1,0],[0,0,1]]
    >>> crystal.cell
    array([[1, 0, 0],
           [0, 1, 0],
           [0, 0, 1]])
    >>> hamiltonian.cell
    array([[3.14159265, 0.        , 0.        ],
           [0.        , 3.14159265, 0.        ],
           [0.        , 0.        , 4.71238898]])

Creation of the Hamiltonian
===========================

The Hamiltonian could be created from the scratch or from the |TB2J|_ file.
For the creation on the basis of the |TB2J|_ file see :py:func:`.read_tb2j_model`.
Here we cover the creation from scratch.

An instance of the :py:class:`.SpinHamiltonian` can be created as:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()

Constructor can take two optional parameters: 

* :py:class:`.Crystal` instance, which defines the structure of the Hamiltonian.

.. doctest::

    >>> import radtools as rad
    >>> lattice = rad.lattice_example("TET")
    >>> crystal = rad.Crystal(lattice)
    >>> hamiltonian = rad.SpinHamiltonian(crystal)

* ``notation``: notation of the Hamiltonian. See :ref:`guide_exchange_notation`.

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian(notation="standard")

When the instance of the :py:class:`.SpinHamiltonian` is created,
it is empty. It need to be filled with bonds and exchange parameters.

Adding atoms
============

Addition of atoms to the :py:class:`.SpinHamiltonian` is directly inherited from the
:py:class:`.Crystal`:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> len(hamiltonian.atoms)
    0
    >>> hamiltonian.add_atom("Cr", spin=1.5, position=(1,0,0))
    >>> hamiltonian.add_atom(rad.Atom("Br"))
    >>> len(hamiltonian.atoms)
    2

However, removing an atom from the :py:class:`.SpinHamiltonian` and from the
:py:class:`.Crystal` is different. When atom is removed from the :py:class:`.SpinHamiltonian`
all bonds, which are connected to it are removed as well. 

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> hamiltonian.add_atom("Cr", spin=1.5, position=(1,0,0))
    >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=1), "Cr", "Cr", (1, 0, 0))
    >>> len(hamiltonian)
    1
    >>> len(hamiltonian.atoms)
    1
    >>> hamiltonian.remove_atom("Cr")
    >>> len(hamiltonian)
    0
    >>> len(hamiltonian.atoms)
    0

Adding bonds
============

The bond is added to the Hamiltonian with :py:meth:`.SpinHamiltonian.add_bond` method:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> Cr = rad.Atom("Cr", spin=1.5)
    >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=1), Cr, Cr, (1, 0, 0))
    >>> hamiltonian["Cr", "Cr", (1, 0, 0)].iso
    1.0

.. note::

    The bond is added to the Hamiltonian, but no atom is explicitly added to the crystal.
    You can skip the atom addition if you are using :py:class:`.Atom` instances, 
    not the string literals in the :py:meth:`.SpinHamiltonian.add_bond` method.
    The atom will be silently added to the system.
    When string literal is used, atom object is extracted from the exchange Hamiltonian, thus
    it has to be explicitly added first:

    .. doctest::

        >>> import radtools as rad
        >>> hamiltonian = rad.SpinHamiltonian()
        >>> Cr = rad.Atom("Cr", spin=1.5)
        >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=1), "Cr", "Cr", (1, 0, 0))
        Traceback (most recent call last):
        ...
        ValueError: No match found for name = Cr, index = None
        >>> hamiltonian.add_atom(Cr)
        >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=1), "Cr", "Cr", (1, 0, 0))
        >>> hamiltonian["Cr", "Cr", (1, 0, 0)].iso
        1.0

Equivalent way to  add bond:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> Cr = rad.Atom("Cr", spin=1.5)
    >>> hamiltonian[Cr, Cr, (1, 0, 0)] = rad.ExchangeParameter(iso=1)
    >>> hamiltonian[Cr, "Cr", (2, 0, 0)] = rad.ExchangeParameter(iso=2)
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    1.0
    >>> hamiltonian[Cr, Cr, (2, 0, 0)].iso
    2.0

.. note::

    If the atom is already added to the Hamiltonian, then you can use both atom`s name and 
    the instance of the atom class. 
    If you have more then one atom with the same name in the Hamiltonian, 
    then you have to use :py:meth:`.Atom.fullname` instead of the name:

    .. doctest::

        >>> import radtools as rad
        >>> hamiltonian = rad.SpinHamiltonian()
        >>> Cr1 = rad.Atom("Cr", position=(0.5, 0, 0))
        >>> Cr2 = rad.Atom("Cr", position=(1, 0, 0))
        >>> hamiltonian.add_atom(Cr1)
        >>> hamiltonian.add_atom(Cr2)
        >>> hamiltonian[Cr1, Cr2, (1, 0, 0)] = rad.ExchangeParameter(iso=2)
        >>> hamiltonian[Cr1, Cr2, (1, 0, 0)].iso
        2.0
        >>> hamiltonian["Cr", "Cr", (1, 0, 0)] = rad.ExchangeParameter(iso=1)
        Traceback (most recent call last):
        ...
        ValueError: Multiple matches found for name = Cr, index = None
        >>> hamiltonian["Cr__1", "Cr__2", (1, 0, 0)] = rad.ExchangeParameter(iso=1)
        >>> hamiltonian["Cr__1", "Cr__2", (1, 0, 0)].iso
        1.0

To remove the bond use :py:meth:`.SpinHamiltonian.remove_bond` method:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> Cr = rad.Atom("Cr", spin=1.5)
    >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=1), Cr, Cr, (1, 0, 0))
    >>> (Cr, Cr, (1, 0, 0)) in hamiltonian
    True
    >>> hamiltonian.remove_bond("Cr", Cr, (1, 0, 0))
    >>> (Cr, Cr, (1, 0, 0)) in hamiltonian
    False

.. hint::

    Note how we used string literal to remove the bond.

Equivalent way to delete bond:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> Cr = rad.Atom("Cr", spin=1.5)
    >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=1), Cr, Cr, (1, 0, 0))
    >>> (Cr, Cr, (1, 0, 0)) in hamiltonian
    True
    >>> del hamiltonian[Cr, Cr, (1, 0, 0)]
    >>> (Cr, Cr, (1, 0, 0)) in hamiltonian
    False

Structure of the Hamiltonian
============================

Main idea of the exchange Hamiltonian structure could be expressed as 

.. code-block:: text

    (atom1, atom2, (i,j,k)) -> exchange_parameter

where ``atom1`` and ``atom2`` are instances of :py:class:`.Atom`, ``(i,j,k)`` is a tuple of
integers, which defines the unit cell, and ``exchange_parameter`` is an instance of
:py:class:`.ExchangeParameter`. 

Exchange Hamiltonian is iterable over bonds:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> Cr = rad.Atom("Cr", spin=1.5)
    >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=1), Cr, Cr, (1, 0, 0))
    >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=2), Cr, Cr, (0, 1, 0))
    >>> for atom1, atom2, (i,j,k), J in hamiltonian:
    ...     print(atom1.fullname, atom2.fullname, (i,j,k), J.iso)
    ... 
    Cr__1 Cr__1 (1, 0, 0) 1.0
    Cr__1 Cr__1 (0, 1, 0) 2.0

It could be check for the presence of the bond:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> Cr = rad.Atom("Cr", spin=1.5)
    >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=1), Cr, Cr, (1, 0, 0))
    >>> (Cr, Cr, (1, 0, 0)) in hamiltonian
    True
    >>> (Cr, Cr, (0, 1, 0)) in hamiltonian
    False

The exchange parameter could be accessed (and set) directly:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> Cr = rad.Atom("Cr", spin=1.5)
    >>> hamiltonian[Cr, Cr, (1, 0, 0)] = rad.ExchangeParameter(iso=1)
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    1.0

``len()`` of the Hamiltonian returns the number of bonds:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> Cr = rad.Atom("Cr", spin=1.5)
    >>> len(hamiltonian)
    0
    >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=1), Cr, Cr, (1, 0, 0))
    >>> len(hamiltonian)
    1
    >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=2), Cr, Cr, (0, 1, 0))
    >>> len(hamiltonian)
    2

Atom as a part of the key
=========================

Atom object or string literal could be used to access atoms in the Hamiltonian. 
String literal is the :py:attr:`.Atom.name` property of the atom if there is only 
one atom with that name in the :py:attr:`SpinHamiltonian.crystal`. Otherwise 
:py:attr:`.Atom.fullname` property of the atom has to be used.

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> hamiltonian.add_atom("Cr", spin=1.5, position=(0,0,0))
    >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=1), "Cr", "Cr", (0, 1, 0))
    >>> for atom1, atom2, (i,j,k), J in hamiltonian:
    ...     print(atom1.fullname, atom2.fullname, (i,j,k), J.iso)
    ...
    Cr__1 Cr__1 (0, 1, 0) 1.0
    >>> hamiltonian.add_atom("Cr", spin=1.5, position=(1,0,0))
    >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=1), "Cr__1", "Cr__1", (0, 1, 0))
    >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=2), "Cr__1", "Cr__2", (0, 0, 0))
    >>> for atom1, atom2, (i,j,k), J in hamiltonian:
    ...     print(atom1.fullname, atom2.fullname, (i,j,k), J.iso)
    ...
    Cr__1 Cr__1 (0, 1, 0) 1.0
    Cr__1 Cr__2 (0, 0, 0) 2.0

If you want to get an :py:class:`.Atom` instance from the :py:class:`.SpinHamiltonian`
you can use :py:meth:`.SpinHamiltonian.get_atom` 
(:py:meth:`.Crystal.get_atom`) method:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> hamiltonian.add_atom(rad.Atom("Cr", spin=1.5, position=(0,0,0)))
    >>> hamiltonian.add_atom(rad.Atom("Cr", spin=1.5, position=(1,0,0)))
    >>> hamiltonian.get_atom("Cr__1").fullname
    'Cr__1'
    >>> hamiltonian.get_atom("Cr__2").fullname
    'Cr__2'

Magnetic atoms
==============

Atoms, which has at least one bond attached to them are called magnetic atoms.
They could be accessed with :py:attr:`.SpinHamiltonian.magnetic_atoms` property:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> hamiltonian.add_atom("Cr", spin=1.5, position=(0,0,0))
    >>> hamiltonian.add_atom("Cr", spin=1.5, position=(1,0,0))
    >>> hamiltonian.magnetic_atoms
    []
    >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=1), "Cr__1", "Cr__1", (0, 1, 0))
    >>> for atom in hamiltonian.magnetic_atoms:
    ...     print(atom.fullname)
    ...
    Cr__1
    >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=1), "Cr__1", "Cr__2", (0, 0, 0))
    >>> for atom in hamiltonian.magnetic_atoms:
    ...     print(atom.fullname)
    ...
    Cr__1
    Cr__2

Filtering the Hamiltonian
=========================

Exchange Hamiltonian could be filtered by distance, template or set of (i,j,k) tuples 
(R vectors).

Use :py:meth:`.SpinHamiltonian.filter` to filter the instance of the
:py:class:`.SpinHamiltonian` and :py:meth:`.SpinHamiltonian.filtered` to get
the filtered copy of the :py:class:`.SpinHamiltonian`:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> hamiltonian.add_atom("Cr", position=(0.25, 0.25, 0))
    >>> hamiltonian.add_atom("Cr", position=(0.75, 0.75, 0))
    >>> bonds = [
    ...     (12, "Cr__1", "Cr__2", (0, 0, 0)),
    ...     (12, "Cr__2", "Cr__1", (0, 0, 0)),
    ...     (12, "Cr__1", "Cr__1", (1, 0, 0)),
    ...     (12, "Cr__1", "Cr__1", (-1, 0, 0)),
    ...     (12, "Cr__2", "Cr__2", (1, 0, 0)),
    ...     (12, "Cr__2", "Cr__2", (-1, 0, 0)),
    ...     (12, "Cr__1", "Cr__1", (0, 2, 0)),
    ...     (12, "Cr__1", "Cr__1", (0, -2, 0)),
    ...     (12, "Cr__2", "Cr__2", (0, 2, 0)),
    ...     (12, "Cr__2", "Cr__2", (0, -2, 0)),
    ...     (12, "Cr__2", "Cr__1", (2, 2, 0)),
    ...     (12, "Cr__1", "Cr__2", (-2, -2, 0)),
    ... ]
    >>> for J, atom1, atom2, R in bonds:
    ...     hamiltonian.add_bond(rad.ExchangeParameter(iso=J), atom1, atom2, R)
    >>> for atom1, atom2, R, J in hamiltonian:
    ...     print(atom1.fullname, atom2.fullname, R, J.iso)
    ...
    Cr__1 Cr__2 (0, 0, 0) 12.0
    Cr__2 Cr__1 (0, 0, 0) 12.0
    Cr__1 Cr__1 (1, 0, 0) 12.0
    Cr__1 Cr__1 (-1, 0, 0) 12.0
    Cr__2 Cr__2 (1, 0, 0) 12.0
    Cr__2 Cr__2 (-1, 0, 0) 12.0
    Cr__1 Cr__1 (0, 2, 0) 12.0
    Cr__1 Cr__1 (0, -2, 0) 12.0
    Cr__2 Cr__2 (0, 2, 0) 12.0
    Cr__2 Cr__2 (0, -2, 0) 12.0
    Cr__2 Cr__1 (2, 2, 0) 12.0
    Cr__1 Cr__2 (-2, -2, 0) 12.0
    >>> filtered_hamiltonian = hamiltonian.filtered(max_distance=1)
    >>> for atom1, atom2, R, J in filtered_hamiltonian:
    ...     print(atom1.fullname, atom2.fullname, R, J.iso)
    ...
    Cr__1 Cr__2 (0, 0, 0) 12.0
    Cr__2 Cr__1 (0, 0, 0) 12.0
    Cr__1 Cr__1 (1, 0, 0) 12.0
    Cr__1 Cr__1 (-1, 0, 0) 12.0
    Cr__2 Cr__2 (1, 0, 0) 12.0
    Cr__2 Cr__2 (-1, 0, 0) 12.0
    >>> filtered_hamiltonian = hamiltonian.filtered(min_distance=2.1)
    >>> for atom1, atom2, R, J in filtered_hamiltonian:
    ...     print(atom1.fullname, atom2.fullname, R, J.iso)
    ...
    Cr__2 Cr__1 (2, 2, 0) 12.0    
    Cr__1 Cr__2 (-2, -2, 0) 12.0
    >>> filtered_hamiltonian = hamiltonian.filtered(R_vector=[(0, 0, 0), (1, 0, 0)])
    >>> for atom1, atom2, R, J in filtered_hamiltonian:
    ...     print(atom1.fullname, atom2.fullname, R, J.iso)
    ...
    Cr__1 Cr__2 (0, 0, 0) 12.0
    Cr__2 Cr__1 (0, 0, 0) 12.0
    Cr__1 Cr__1 (1, 0, 0) 12.0
    Cr__2 Cr__2 (1, 0, 0) 12.0
    >>> filtered_hamiltonian = hamiltonian.filtered(template=[("Cr__1", "Cr__2", (0, 0, 0))])
    >>> for atom1, atom2, R, J in filtered_hamiltonian:
    ...     print(atom1.fullname, atom2.fullname, R, J.iso)
    ...
    Cr__1 Cr__2 (0, 0, 0) 12.0

.. hint::

    Filtering options may be combined together. Only the bonds, which satisfy all the
    conditions will be included in the filtered Hamiltonian.


Energy
======

Ferromagnetic energy of the Hamiltonian could be calculated with
:py:meth:`.SpinHamiltonian.ferromagnetic_energy` method:

.. note::

    The notation of the Hamiltonian has to be defined in order to calculate the energy.

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.SpinHamiltonian()
    >>> hamiltonian.add_atom("Cr", spin=1.5, position=(0,0,0))
    >>> hamiltonian.add_atom("Cr", spin=1.5, position=(1,0,0))
    >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=1), "Cr__1", "Cr__1", (0, 1, 0))
    >>> hamiltonian.add_bond(rad.ExchangeParameter(iso=2), "Cr__1", "Cr__2", (0, 0, 0))
    >>> hamiltonian.ferromagnetic_energy()
    Traceback (most recent call last):
    ...
    radtools.exceptions.NotationError:
    ...
    >>> hamiltonian.notation = "standard"
    >>> hamiltonian.ferromagnetic_energy()
    -13.5

Saving the Hamiltonian
======================

The Hamiltonian could be saved in as a text file with 
:py:meth:`.SpinHamiltonian.dump_txt` method. 

It can be serialized with :py:meth:`.SpinHamiltonian.dump_pickle` method.
