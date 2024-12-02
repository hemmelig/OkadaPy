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
import os
import glob
import sys
from unittest.mock import Mock

import sphinx_rtd_theme
from importlib.metadata import version as _version


sys.path.insert(0, os.path.abspath('..'))

MOCK_MODULES = ["okada.core.libnames"]
for mod_name in MOCK_MODULES:
    sys.modules[mod_name] = Mock()

# -- Project information -----------------------------------------------------
project = "OkadaPy"
copyright = "2023, Conor A. Bacon"
author = "Conor A. Bacon"

# The full version, including alpha/beta/rc tags
release = _version("okada")

# -- General configuration ---------------------------------------------------
# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx_rtd_theme",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosectionlabel",
]

# Set master doc
master_doc = "index"

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# generate automatically stubs
autosummary_generate = glob.glob("submodules" + os.sep + "*.rst")

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_logo = 'img/OkadaPyLogo.png'
html_static_path = []

# Make version number accessible within individual rst files
rst_epilog = """
.. |Version| replace:: {release}
""".format(release=release)
