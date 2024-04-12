.. _contribute_branches:

************
Branch rules
************

Branches of the upstream repository
===================================

In the upstream repository the releases are made from the ``stable`` branch. The
``stable`` branch is the default branch of the upstream repository.

For each release we create a tag named ``vMajor.Minor.Rest`` in the upstream repository.

Starting from the version 0.5, once new ``Minor`` release is made, the ``Major.Minor``
branch is created for the latest release of the previous ``Minor`` version. For example,
when the ``v0.6.0`` is released, the ``0.5`` branch is created from the ``v0.5.29`` tag.
(Provided that the previous release was ``v0.5.29``.)

The ``Major.Minor`` branches are used to backport bug fixes to the previous releases.

Branches of your local repository
=================================

When you clone the upstream repository, you will have the ``stable`` branch checked out by default.

When you want to contribute to the project, you should create a new branch from the ``stable``
branch and name it according to the feature you are working on. If you are fixing a bug, it
is wise to name the branch ``bugfix-issue_number``.

Branches of origin repository
=============================

At any moment of your work on the feature or bug fix, you can push the changes to your
forked repository. For the first push, you should create a new branch in your forked
repository. For example if you named your local branch as ``new-cool-feature``, you can type::

    git push -u origin new-cool-feature

It will create a new branch in your forked repository and link your local branch to it.
For the next pushes, you can simply type ``git push``.
