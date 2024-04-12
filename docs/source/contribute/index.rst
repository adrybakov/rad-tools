.. _contribute:

*******************
Contributor`s guide
*******************


We welcome the contribution to the package!

If you're interested in seeing who has already contributed to this project,
please visit our :ref:`Contributors page <contribute_contributors>`.
We appreciate all contributions and look forward to see your name on that list!

It is not necessary to be a programmer to contribute.
You can help us with the documentation, :ref:`new features <contribute_feature>`
and :ref:`finding bugs <contribute_bug>`.


Contribution to the source code is summarized below.
We assume that you have an account on `<https://github.com>`_
and familiar with `Git <https://git-scm.com/>`_.

Development workflow
====================

For the detailed explanation of the development workflow see the table of content at the
end of the page or click on the links in the short summary below.

Fork and clone
--------------

* Go to the |RAD-repo|_ (:ref:`upstream <contribute_origin-upstream_upstream>`) and
  click on the "Fork" button. Now you have your own copy of the RAD-tools repository
  under your GitHub account (:ref:`origin <contribute_origin-upstream_origin>`).
* Clone your fork to your local machine:

  - If you are using ssh-key::

      git clone git@github.com:your-username/rad-tools.git

  - If you are not using ssh-key::

      git clone https://github.com/your-username/rad-tools.git

* Change the directory::

    cd rad-tools

  Now you are in the root folder of you local repository
  (:ref:`local <contribute_origin-upstream_local>`)

* Add the :ref:`upstream <contribute_origin-upstream>` repository to your
  :ref:`local <contribute_origin-upstream_local>` repository::

    git remote add upstream https://github.com/adrybakov/rad-tools.git

* Pull the latest changes from the RAD-tools repository if necessary::

    git pull upstream stable

  or::

    git fetch upstream
    git merge upstream/stable

Set up the environment
----------------------

We recommend to use virtual environment. Once the virtual environment is created,
you can install requirements:

* Package dependencies::

    pip install -r requirements.txt

* For the package development::

    pip install -r requirements-dev.txt

* For the documentation::

    pip install -r docs/requirements.txt

* For the tests::

    pip install -r utest/requirements.txt

.. note::
  For the linux and OSX systems there is a scenario defined.
  It installs all requirements. Note: it does NOT create a virtual environment::

    make requirements

Enable pre-commit
-----------------

We use `pre-commit <https://pre-commit.com/>`_ to enforce some rules on the code style
before each commit.
To enable it, run the following command::

  pre-commit install

Now, every time you commit the code, pre-commit will check it for you. If some checks fail,
pre-commit will automatically fix them and abort the commit. You need to add the fixed files
to the staging area and commit again.

.. hint::
  If you want to run pre-commit manually, you can use the following command::

    pre-commit run --all-files

Develop your contribution
-------------------------

* Create a :ref:`dedicated branch <contribute_branches>` for your feature,
  that you are going to develop::

    git checkout -b feature-name

* Develop your contribution. Commit your progress locally
  (`git-add <https://git-scm.com/docs/git-add>`_
  and `git-commit <https://git-scm.com/docs/git-commit>`_).
  Use |good-commit-messages|_. Write :ref:`tests <contribute_tests>`.
  Write :ref:`documentation <contribute_documentation>`.

Submit your contribution
------------------------

* Push the changes to your forked repository::

    git push origin feature-name

* Go to your forked repository on GitHub and click on the
  green "Compare & pull request" button.
  Describe your contribution and submit the pull request.
  Please mention the issue number if it is related to any.

Review and merge
----------------

* Once the pull request is submitted, the code will be reviewed.
  If there are any comments, please fix them. You can make the changes locally, commit
  them and push to the same branch of :ref:`origin <contribute_origin-upstream_origin>`
  repository and they will be added to the pull request automatically.
* Once the pull request is approved, it will be merged to the
  `stable <https://github.com/adrybakov/rad-tools>`_ branch.



Development process in details
==============================

.. toctree::
   :hidden:

   contributors

.. toctree::
   :maxdepth: 2

   origin-upstream
   branches
   features
   bugs
   documentation
   tests
   script
