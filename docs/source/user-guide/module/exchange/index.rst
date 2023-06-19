.. _rad-tools_exchange:

********************
Exchange Hamiltonian
********************

For the full reference see :ref:`api_exchange`.

.. currentmodule:: radtools


:py:class:`.ExchangeHamiltonian` is build on some :py:class:`.Crystal`, 
which defines the structure (lattice + atoms). It stores bonds between :py:class:`.Atom` 
:math:`i`, which is located in :math:`(0, 0, 0)` unit cell and :py:class:`.Atom` :math:`j`,
which is located in :math:`(i, j, k)` unit cell and corresponding exchange parameters.

The notation of the exchange Hamiltonian is an important issue, since without clear notation 
exchange parameters does not make much sense. :py:class:`.ExchangeHamiltonian` can support any 
notation out of the most common ones. Detailed description of the notation is given in 
correspondent section: :ref:`rad-tools_exchange-notation`.

The main building block of the exchange Hamiltonian is exchange parameter 
(:math:`\boldsymbol{J}`). It is separated in a separate class :py:class:`.ExchangeParameter`.
For the guide on exchange parameter see:

.. toctree::
    :maxdepth: 1
    
    parameter



.. _rad-tools_exchange-notation:

Notation
========

Here is the Hamiltonian in the notation, which is considered to be the standard for the
RAD-tools:

.. math::

    H = -\sum_{i,j} \hat{\boldsymbol{S}}_i \cdot \boldsymbol{J}_{i,j} \hat{\boldsymbol{S}}_j

It has double-counting, spins are not normalized, and the exchange parameter is
positive for ferromagnetic order.

.. note::
    The "standard" here does not mean that the :py:class:`.ExchangeHamiltonian` is
    always has this notation. For example when :py:class:`.ExchangeHamiltonian` is
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
exchange Hamiltonian. During the creation of the ExchangeHamiltonian object the notation 
is deliberately not defined, because it depends on your interpretation. Therefore,
the notation has to be defined explicitly by you. If the notation is not defined, 
and the you are trying to use the properties and methods, which expect the notation
to be defined, then the :py:exc:`.NotationError` is raised.

The notation could be defined in two ways:

* By setting each individual property. For example:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.ExchangeHamiltonian()
    >>> hamiltonian.double_counting = True
    >>> hamiltonian.spin_normalized = False
    >>> hamiltonian.factor_one_half = False
    >>> hamiltonian.factor_two = True
    >>> hamiltonian.minus_sign = True

* By setting the :py:attr:`.notation` property directly. For example:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.ExchangeHamiltonian()
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
parameters. See :ref:`Examples <rad-tools_exchange-notation-examples>`.

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
    >>> hamiltonian = rad.ExchangeHamiltonian()
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

.. _rad-tools_exchange-notation-examples:

Examples
--------
Setting one of the predefined notations:

.. doctest::

    >>> import radtools as rad
    >>> hamiltonian = rad.ExchangeHamiltonian()
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
    >>> hamiltonian = rad.ExchangeHamiltonian()
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
    >>> hamiltonian = rad.ExchangeHamiltonian()
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
    >>> hamiltonian = rad.ExchangeHamiltonian()
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

