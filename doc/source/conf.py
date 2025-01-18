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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

import textwrap

# -- Project information -----------------------------------------------------

project = 'Idefix'
copyright = 'Geoffroy Lesur et al.'
author = 'Geoffroy Lesur'

# The full version, including alpha/beta/rc tags
release = '2.2.00'



# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx_rtd_theme",
    'sphinx_git',
    "breathe",
    "exhale",
    "m2r2",
    "sphinx_copybutton"
    ]

source_suffix = [".rst", ".md"]

cpp_id_attributes=["KOKKOS_FORCEINLINE_FUNCTION","KOKKOS_INLINE_FUNCTION","KOKKOS_RESTRICT"]

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
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

def setup(app):
    app.add_css_file('my_theme.css')

##############
# Breathe/Exhale options

# Setup the breathe extension
breathe_projects = {
    "Idefix": "./xml"
}
breathe_default_project = "Idefix"

# Setup the exhale extension
exhale_args = {
    # These arguments are required
    "containmentFolder":     "./api",
    "rootFileName":          "library_root.rst",
    "rootFileTitle":         "Idefix API",
    # Suggested optional arguments
    "createTreeView":        False,
    # TIP: if using the sphinx-bootstrap-theme, you need
    # "treeViewIsBootstrap": True,
    "doxygenStripFromPath": "../../src",
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin":    textwrap.dedent('''
            INPUT = ../../src
            EXCLUDE = ../../src/kokkos
            ''')
}

# Tell sphinx what the primary language being documented is.
primary_domain = 'cpp'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'
