# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import subprocess
import sys
import re
# sys.path.insert(0, os.path.abspath('.'))

# -- Project information -----------------------------------------------------

project = u'fit2x'
copyright = u'2019, Thomas Peulen & Seidel Lab (Heinrich-Heine University)'
author = u'Thomas-Otavio Peulen'

# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
needs_sphinx = '1.8.5'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'matplotlib.sphinxext.plot_directive',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.autodoc',
    'sphinxcontrib.bibtex',
    'recommonmark',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.intersphinx',
    'nbsphinx',
    'sphinx_gallery.gen_gallery'
]


class SubSectionTitleOrder:
    """Sort example gallery by title of subsection.
    Assumes README.rst exists for all subsections and uses the subsection with
    dashes, '---', as the adornment.
    """
    def __init__(self, src_dir):
        self.src_dir = src_dir
        self.regex = re.compile(r"^([\w ]+)\n-", re.MULTILINE)

    def __repr__(self):
        return '<%s>' % (self.__class__.__name__,)

    def __call__(self, directory):
        src_path = os.path.normpath(os.path.join(self.src_dir, directory))

        # Forces Release Highlights to the top
        if os.path.basename(src_path) == "release_highlights":
            return "0"

        readme = os.path.join(src_path, "README.rst")

        try:
            with open(readme, 'r') as f:
                content = f.read()
        except FileNotFoundError:
            return directory

        title_match = self.regex.search(content)
        if title_match is not None:
            return title_match.group(1)
        return directory


# This is largely from
# https://github.com/scikit-learn/scikit-learn/blob/master/doc/conf.py
sphinx_gallery_conf = {
    'doc_module': 'fit2x',
    'backreferences_dir': os.path.join('modules', 'generated'),
    'show_memory': False,
    'examples_dirs': ['../examples'],
    'gallery_dirs': ['auto_examples'],
    'subsection_order': SubSectionTitleOrder('../examples'),
    'inspect_global_variables': False,
    'remove_config_comments': True,
    'filename_pattern': '/plot_',
    'ignore_pattern': r'__init__\.py',
}

nbsphinx_execute_arguments = [
    "--InlineBackend.figure_formats={'svg', 'pdf'}",
    "--InlineBackend.rc={'figure.dpi': 96}",
]

# numpydoc
numpydoc_attributes_as_param_list = False

# matplotlib plot directive
plot_include_source = False
plot_formats = [("png", 90)]
plot_html_show_formats = True
plot_html_show_source_link = True
plot_pre_code = """import numpy as np
import fit2x
import pylab as p
"""


breathe_projects = {}
# latex and breathe do not play very well together. Therefore
# breathe is only used for the webpage.
# compatible with readthedocs online builder and local builder
if sys.argv[0].endswith('sphinx-build') and \
        ('html' in sys.argv or sys.argv[-1] == '_build/html'):
    subprocess.call('doxygen', shell=True)
    breathe_projects['fit2x'] = './_build/xml'
    breathe_default_project = "fit2x"
    extensions += ['breathe']


autosectionlabel_prefix_document = True
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#

#source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [u'_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = None


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

# Read the docs style:
on_rtd = os.environ.get('READTHEDOCS') == 'True'
if not on_rtd:
    import sphinx_bootstrap_theme
    # Activate the theme.
    html_theme = 'bootstrap'
    html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()
else:
    try:
        import sphinx_rtd_theme
    except ImportError:
        pass  # assume we have sphinx >= 1.3
    else:
        html_theme_path = sphinx_rtd_theme.get_html_theme_path()
    html_theme = 'sphinx_rtd_theme'


def setup(app):
    app.add_css_file("fix_rtd.css")


# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'fit2xdoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'fit2x.tex', u'fit2x Documentation', u'Thomas-Otavio Peulen', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).ƒre
man_pages = [
    (master_doc, 'fit2x', u'fit2x Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'fit2x', u'fit2x Documentation',
     author, 'fit2x', 'One line description of project.',
     'Miscellaneous'),
]


# -- Options for Epub output -------------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#
# epub_identifier = ''

# A unique identification for the text.
#
# epub_uid = ''

# A list of files that should not be packed into the epub file.
epub_exclude_files = ['search.html']


# -- Extension configuration -------------------------------------------------
