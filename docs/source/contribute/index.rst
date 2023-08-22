.. _contribute:

*******************
Contributor`s guide
*******************

We welcome the contribution to the package! 

It is not necessary to be a programmer to contribute. 
You can help us with the documentation, :ref:`new features <contribute_feature>` 
and :ref:`finding bugs <contribute_bug>`.

.. toctree::
   :hidden:

   features
   bugs
   documentation

Contribution to the source code is summarized below.

Fork
====
The first step is to fork the repository on GitHub. 
Go to |RAD-repo|_ and click on the "Fork" button.

Clone
=====
Clone you forked repository to your local machine.

If you are using ssh-key:

.. code-block:: bash

   git clone https://github.com/your-username/rad-tools.git

If you are not using ssh-key:

.. code-block:: bash

   git clone git@github.com:your-username/rad-tools.git

Change the directory:

.. code-block:: bash

   cd rad-tools

Environment
===========

Set up the necessary environment. We recommend to use virtual environment.
Once the virtual environment is created, you can install requirements:

* For the package development:

.. code-block:: bash

   pip install -r requirements.txt

* For the docs:

.. code-block:: bash

   pip install -r docs/requirements.txt

* For the tests:

.. code-block:: bash

   pip install -r utest/requirements.txt

Develop
=======

Develop your contribution. Commit your progress locally. Use meaningful commit messages.
If you want to create new script for command line interface, then check out :ref:`contribute_script` guide.

.. toctree::
   :hidden:

   script

Test
====

Write unit tests for the developed code. We recommend to write tests before the code.
Test files are located in the "utest" directory. the structure of the source code directory 
("radtools") and the test directory ("utest") is loosely the same.

Run the tests locally:

.. code-block:: bash

   make test

Docs
====

If you have developed a new feature, please add the description to the docs. 
Make sure to follow :ref:`docs_guide`.

Build the docs locally:

.. code-block:: bash

   make html

Doctest
=======

Run doctest locally:

.. code-block:: bash

   make doctest

Pull request
============

Once you are ready with your contribution, 
push your changes to your forked repository and create a pull request.
Please add the description to the pull request. 
Mention the issue number if it is related to any.
