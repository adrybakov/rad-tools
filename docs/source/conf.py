from os.path import abspath
import sys
from rad_tools import __version__
sys.path.insert(0, abspath('..'))


# Project information
project = 'rad-tools'
copyright = '2022, Andrey Rybakov'
author = 'Andrey Rybakov'

# Project version
release = __version__


# -- General configuration
extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.viewcode',
    'sphinx.ext.extlinks',
    'numpydoc',
    'sphinx_copybutton',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.mathjax',
    'sphinx_design',
    'sphinx.ext.doctest',
]

autosummary_generate = True

smartquotes = False

# Extlinks
extlinks = {'ANSI': (
    'https://en.wikipedia.org/wiki/ANSI_escape_code%s', None),
    'examples': (
    'https://github.com/adrybakov/rad-tools/tree/stable/docs/examples/%s', '%s'),
    'DOI': (
    'https://doi.org/%s', 'DOI: %s'), }

# todo
todo_include_todos = True

# Solve the problem with warnings
numpydoc_class_members_toctree = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

#pydata bugfix
icon_links = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = 'pydata_sphinx_theme'
# html_theme = 'alabaster'
# html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']
html_css_files = ["rad-tools.css"]

# Logo
# html_logo = 'img/logo-black.png'
html_title = f'{project} {release}'
html_favicon = 'img/favicon.png'
html_theme_options = {
    "github_url": "https://github.com/adrybakov/rad-tools",
    "twitter_url": "https://twitter.com/adrybakov",
    "show_nav_level": 2,
    "use_edit_page_button": True,
    "navbar_end": ["theme-switcher.html", "navbar-icon-links.html"],
    "logo": {
        "image_light": "logo_black.png",
        "image_dark": "logo_white.png",
    },
}

html_context = {
    "default_mode": "light",
    "display_github": True,  # Integrate GitHub
    "github_user": "adrybakov",  # Username
    "github_repo": "rad-tools",  # Repo name
    "github_version": "master",  # Version
    "doc_path": "docs/source",  # Path in the checkout to the docs root
    'docsearch_disabled': False,
}


# Custom variables with access from .rst files
variables_to_export = [
    "project",
    "copyright",
    "release",
]
frozen_locals = dict(locals())
rst_epilog = '\n'.join(
    map(lambda x: f".. |{x}| replace:: {frozen_locals[x]}", variables_to_export))
del frozen_locals
