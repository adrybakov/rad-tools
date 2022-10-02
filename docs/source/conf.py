# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
from os.path import abspath
import sys
from rad_tools import __version__
sys.path.insert(0, abspath('..'))


# -- Project information -----------------------------------------------------

project = 'rad-tools'
copyright = '2022, Andrey Rybakov'
author = 'Andrey Rybakov'

# The full version, including alpha/beta/rc tags
release = __version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'numpydoc',
    'sphinx_rtd_theme',
    'sphinx_copybutton',
]

# Solve the problem with warnings
numpydoc_class_members_toctree = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

# html_theme = 'sphinx_book_theme'
# html_theme = 'alabaster'
html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']

# Logo
html_logo = 'logo.png'
html_theme_options = {
    'logo_only': True,
    'display_version': False,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'style_nav_header_background': '#4E7F6C',
    # Toc options
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False}

# Integration of GitHub
html_context = {
    "display_github": True, # Integrate GitHub
    "github_user": "adrybakov", # Username
    "github_repo": "rad-tools", # Repo name
    "github_version": "master", # Version
    "conf_py_path": "/docs/", # Path in the checkout to the docs root
}

# Variables with access from .rst files
variables_to_export = [
    "project",
    "copyright",
    "release",
]
frozen_locals = dict(locals())
rst_epilog = '\n'.join(map(lambda x: f".. |{x}| replace:: {frozen_locals[x]}", variables_to_export))
del frozen_locals

