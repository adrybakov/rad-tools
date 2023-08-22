.. _contribute_script:

**********
New script
**********

Each script is meant to be an interface to the Python library. It is typically implemented 
as a small ``manager()`` function in the :ref:`api_score` module.

Workflow
========

* Create templates for script and docs:

.. code-block:: bash

    make new-script NAME=script-name

* Implement the manager() function.

* Write docs for the script.

* Run post-processors, which creates Arguments section in docs and create_parser() function in script.

.. code-block:: bash

    make generate-script-docs SCRIPT=script-name

Documentation for the manager function
======================================

.. note::
    For correct post-processing your python files have to be formatted with |black|_.

Check the documentation of other scripts for examples.

You need to follow guidelines of |numpydoc|_ and some additional rules:

* Each parameter has to have "Console argument:" in the description.
* Description goes to the help message of the parser.
* Description has to be one line long.
* Description and long description goes to the documentation.
* Console arguments have to be separated by "/"
* "list of str" means that nargs="*", type="str" passed to ``parser.add_argument()``.
* "list of 2 str" means that nargs=2, type="str" passed to ``parser.add_argument()``.
* If default value is ``None``, then word "optional" word appears in the documentation.
* If default value is not ``None``, then word "default: Value" word appears in the documentation and ``default=Value`` passed to ``parser.add_argument()``.
* If argument is positional, then ``required=True`` passed to ``parser.add_argument()``.


.. code-block:: python

    def manager(input_file, optional_param_1=None, optional_param_2=None, default_param=1):
        r"""
        ...

        Parameters
        ==========
        input_file : str
            Description.

            Long description.

            Console argument: ``-if`` / ``--input-file``
        optional_param_1 : list of str, optional
            Description.

            Long description.

            Console argument: ``-op1`` / ``--optional-param-1``
        optional_param_2 : list of 2 str, optional
            Description.

            Long description.

            Console argument: ``-op2`` 
        default_param : int, default 34
            Description.

            Long description.

            Console argument: ``-dp`` / ``--default-param``
        """

