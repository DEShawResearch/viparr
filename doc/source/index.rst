
Viparr Documentation
====================

Viparr is a tool that applies forcefield parameters to a chemical system
file, producing an output file that can be used to perform molecular
dynamics simulations with Anton or Desmond.

It contains a collection of smaller command-line tools to handle a number
of forcefield-related tasks, as well as a Python library enabling easy
customization and extension of certain aspects of viparr's execution.

*New in viparr/4.7.0*: The way forcefield directories are discovered
has changed.  Previously, forcefields specified by name (rather than
by path) needed to be located in a single directory given by one of the
environment variables `VIPARR3_FFDIR` or `VIPARR4_FFDIR`.  Forcefields are
now searched using the environment variable `VIPARR_FFPATH`, which is
a colon-delimited list of directories.   Use `viparr-ff/2.1.0` or later
with this later later versions of viparr.

.. toctree::
   :maxdepth: 2

   viparr_normal
   tools
   python_api
   release_notes.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

