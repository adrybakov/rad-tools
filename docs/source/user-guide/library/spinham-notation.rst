.. _library_spinham_notation:

****************************
Notation of spin Hamiltonian
****************************

The notation, which is considered to be the standard for the RAD-tools:

.. math::

    H = -\sum_{i,j} \hat{\boldsymbol{S}}_i \cdot \boldsymbol{J}_{i,j} \hat{\boldsymbol{S}}_j

It has double-counting, spins are not normalized, and the exchange parameter is
positive for ferromagnetic order.

.. note::
    The "standard" here does not mean that the :py:class:`.SpinHamiltonian`
    always has this notation. For example when :py:class:`.SpinHamiltonian` is
    read from |TB2J|_ file (:py:func:`.read_tb2j_model`) it has the notation of TB2J.

There are five properties, that has to be defined in order to describe the 
spin Hamiltonian`s notation:

* :py:attr:`.double_counting`
    Whether both pairs :math:`(i, j)` and :math:`(j, i)` are included in the sum. 

    If double counting is avoided it is usually indicated under the sum sign:

    .. math::

        \sum_{i < j} \text{ or } \sum_{i > j}

    When double counting is implied there are no indication in the sum: 

    .. math::

        \sum_{ij}

    However, some authors are using 
    
    .. math::
        
        \sum_{<ij>} 
    
    as an indication of the avoided double counting, which we find confusing 
    and discourage you to use it in this way. In textbooks this indication usually 
    implies that only near-neighbors are included in the sum, which is not the same as 
    avoiding double counting (consecutively :math:`\sum_{<<ij>>}` is used to indicate 
    next-nearest neighbors and so on).

    .. note::

        Indication

        .. math::

            \sum_{i \ne j}

        Does not mean that double counting is avoided. 
        It specifies that there is no exchange between the same atom.

    
* :py:attr:`.spin_normalized`
    :math:`\boldsymbol{S}_i` in the spin Hamiltonian is viewed as spin operator 
    or as classical spin vector. It could be normalized to 1 or not. 
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
spin Hamiltonian. During the creation of the :py:class:`.SpinHamiltonian` object the 
notation is deliberately not defined, because it depends on your interpretation. 
Therefore, the notation has to be defined explicitly by you. If the notation is not 
defined, and the you are trying to use the properties and methods, which expect the 
notation to be defined, then the :py:exc:`.NotationError` is raised. 
See :ref:`Examples <library_spinham_notation-examples>` for the usage.


Conversion of the notation
==========================

In the figure below we illustrate the problem of the notation. The picture describe the 
whole notation tree, if the no minus sign is present (:py:attr:`.SpinHamiltonian.minus_sign`
is ``False``). Click on the picture to enlarge it.

.. figure:: ../../img/notation-tree.png
    :target: ../../_images/notation-tree.png

In the picture we assume that there are only two centers :math:`i` and :math:`j`, which are
connected by the exchange parameter :math:`\boldsymbol{J}_{i,j}`. The alignment of the
spins is ferromagnetic. Then the numbers on the right side of the picture are the
energies of that system in the corresponding notation. For each notation the exchange
parameter is set to :math:`1`.

These numbers provide a way for you to convert the exchange parameters from one notation
to another. The conversion is done based on the fact that any notation should give the same energy
for the same system. Note, that in the picture the systems are different since the exchange 
parameters are equal to one in each notation, but the energies are different.

The rule for the conversion from notation one (N1) to the notation two (N2) is the following:

* Follow the notation tree for each notation. Obtain factor one :math:`f_1` and factor two :math:`f_2`.

* If in the notation one minus sign is written (i.e. positive exchange parameter means ferromagnetic alignment), multiply the factor :math:`f_1` by :math:`-1`.

* If in the notation two minus sign is written (i.e. positive exchange parameter means ferromagnetic alignment), multiply the factor :math:`f_2` by :math:`-1`.

* compute conversion factor :math:`f` for the exchange parameters as

.. math::

    f = \dfrac{f_1}{f_2}

* Multiply exchange parameters in the notation one by the conversion factor :math:`f`.

The result are the exchange parameters in the notation two.



Predefined notations
====================

There are three predefined notations in the RAD-tools. Each predefined notation is a 
tuple of five ``bool``, which correspond to the five properties of the notation.

.. hint::
    Order: (double counting, spin normalized, factor 1/2, factor 2, minus sign).

* Standard
    (True, False, False, False, True)

    .. math::
        H = -\sum_{i,j} \hat{\boldsymbol{S}}_i \cdot \boldsymbol{J}_{i,j} \cdot \hat{\boldsymbol{S}}_j

    where double counting is present (:math:`ij` and :math:`ji` is in the sum).
    Spin vectors are **not** normalized.
* |TB2J|_
    (True, True, False, False, True)

    .. math::
        H = -\sum_{i,j} \hat{\boldsymbol{S}}_i \cdot \boldsymbol{J}_{i,j} \cdot \hat{\boldsymbol{S}}_j

    where double counting is present (:math:`ij` and :math:`ji` is in the sum).
    Spin vectors are normalized to 1.
* SpinW
    (True, False, False, False, False)

    .. math::
        H = \sum_{i,j} \hat{\boldsymbol{S}}_i \cdot \boldsymbol{J}_{i,j} \cdot \hat{\boldsymbol{S}}_j

    where double counting is present (:math:`ij` and :math:`ji` is in the sum).
    Spin vectors are **not** normalized.

See :ref:`Examples <library_spinham_notation-examples>` for the usage.