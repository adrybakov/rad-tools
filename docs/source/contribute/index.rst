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

For the detailed explanation of the development workflow, please visit
the corresponding links below.

1. Fork and clone.

   * Go to the |RAD-repo|_ and click on the "Fork" button.
     Now you have your own copy of the RAD-tools repository in your GitHub account.
   * Clone your copy of the repository to your local machine:

     - If you are using ssh-key::

         git clone git@github.com:your-username/rad-tools.git

     - If you are not using ssh-key::

         git clone https://github.com/your-username/rad-tools.git

   * Change the directory::

      cd rad-tools
   * Add the upstream repository::

      git remote add upstream https://github.com/adrybakov/rad-tools.git

2. Set up the environment.

   We recommend to use virtual environment. Once the virtual environment is created,
   you can install requirements:

   * For the package development::

      pip install -r requirements.txt
   * For the docs::

      pip install -r docs/requirements.txt
   * For the tests::

      pip install -r utest/requirements.txt

3. Develop the contribution.

   * Create a dedicated branch for your feature, that you are going to develop::

      git checkout -b feature-name

   * Develop your contribution. Commit your progress locally
     (`git-add <https://git-scm.com/docs/git-add>`_
     and `git-commit <https://git-scm.com/docs/git-commit>`_).
     Use meaningful commit messages. Write `tests <contribute_tests>`.
     Write `documentation <contribute_docs>`.

4. Submit your contribution.

   * Push the changes to your forked repository::

      git push origin feature-name

   * Go to your forked repository on GitHub and click on the
     green "Compare & pull request" button.
     Describe your contribution and submit the pull request.
     Please mention the issue number if it is related to any.

5. Review and merge.

   * Once the pull request is submitted, the code will be reviewed.
     If there are any comments, please fix them.
   * Once the pull request is approved, it will be merged to the
     `stable <https://github.com/adrybakov/rad-tools>`_ or
     `dev <https://github.com/adrybakov/rad-tools/tree/dev>`_ branch.


Development process in details
==============================

.. toctree::
   :hidden:

   contributors

.. toctree::
   :maxdepth: 2

   features
   bugs
   documentation
   tests
   script
