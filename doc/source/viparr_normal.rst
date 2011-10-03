````````````````````
Executing ``viparr``
````````````````````

``viparr`` executable
=====================

::

  $ viparr input_system output.dms [(-f ffname | -d ffdir) [-m ffpatch]* \
      [-a ffpatch]* ]* [-s selection] [--reorder-ids] [--rename-atoms] \
      [--rename-residues] [--without-constraints] [--without-fix-masses] \
      [--non-fatal] [--verbose-plugins]

``viparr`` applies forcefield parameters to a chemical system. This occurs in
four main steps:

#. Any forcefield patches specified with -m or -a are merged with the last
   preceding -d or -f forcefield.

#. The resulting forcefields are used to typify the selected atoms of the
   chemical system. Each bonded molecule in the atom selection is typified by
   the first forcefield that can match all of its atoms.

#. The parameters of each forcefield are applied to the atoms typified by that
   forcefield using the plugins specified in the forcefield's rules file.

#. Post-processing steps are applied to build constraints, equate masses of
   atoms of the same element, and reorder atom IDs. These steps may be
   controlled by the available command-line flags.

The ``input_system`` may be in MAE or DMS format. To use template-based
forcefields with residue templates, atoms in the input system must be correctly
demarcated into residues using a combination of fragment, chain name/ID, and
residue name/ID. Residues in different fragments and different chains can have
the same residue number. ``viparr`` uses atomic number and bond structure (graph
isomorphism) to match residues to forcefield templates; atom, residue, and
template names are ignored in the matching process.

During parameter matching, the plugins determine how the matching is to be done
(i.e. which permutations of an atom tuple can be matched to each row of the
forcefield's parameter table, whether to match bond orders if they are present
in the parameter table, how to apply the matches to the output system, etc.).
An atom tuple will be matched to the parameter row specifying its atom types
exactly or the first row that matches its atom types with at least one wildcard
'*' type.

To avoid ambiguity, the VDW functional form and VDW combine rule, if defined in
multiple forcefields or forcefield patches, must be identical. If ``selection``
does not contain the entire chemical system, then the existing VDW functional
form and VDW combine rule of the system are preserved and must match those
defined in the forcefields as well. If no VDW functional forms or combine rules
are specified, defaults of ``lj12_6_sig_epsilon`` and ``geometric`` are assumed.

Unless ``--reorder-ids`` is specified, ``viparr`` preserves the input atom IDs
in the output DMS file. ``viparr`` also preserves residue numbering and chain
numbering, as well as atom and residue names by default. Existing pseudo atoms
with atomic number 0, bonds to these pseudo atoms, forcefield parameters, and
constraints in the input DMS file corresponding to the given ``selection`` are
discarded, and new pseudo atoms, pseudo bonds, parameters, and constraints are
added for this ``selection`` in ``output.dms``. If ``selection`` does not
contain the entire system, then pseudo atoms, pseudo bonds, parameters, and
constraints for the molecules not in the selection are preserved in
``output.dms``.

Options:

.. program:: viparr

.. cmdoption:: -f ffname

   Forcefield directory found in VIPARR_FFPATH.

.. cmdoption:: -d ffdir        

   Explicit forcefield directory.

.. cmdoption:: -m ffpatch

   Forcefield patch to be merged with the previous forcefield.

.. cmdoption:: -a ffpatch

   Like -m ffpatch, except that elements of ffpatch cannot overwrite those of the original forcefield.

.. cmdoption:: -s selection

   Selection of atoms to parameterize; must correspond to a subset of complete connected molecules. Default is 'all'.

.. cmdoption:: --reorder-ids

   Reorder atom IDs such that IDs of added pseudo particles are next to their parent atoms.
 
.. cmdoption:: --rename-atoms, --rename-residues

   Overwrite atom and residue names in the system with those from the template files.
 
.. cmdoption:: --without-constraints

   Do not build AHn and HOH bond and angle constraints at the end of viparr execution.

.. cmdoption:: --without-fix-masses

   Do not set masses of all atoms of the same element to their median value.

.. cmdoption:: --non-fatal

   Do not exit in error if a required term cannot be matched by any parameter in the forcefield.

.. cmdoption:: --verbose-plugins

   Print a detailed list of which external plugin libraries were loaded and ignored.

Forcefield format
=================

A ``viparr`` forcefield is a directory of files specified in JSON
format. (See http://json.org for details.) The forcefield format for
``viparr`` versions 3.0 and later differs from that for ``viparr``
version 1.*.

Loading the ``viparr-ff/*/data`` module will append to the environment variable
VIPARR_FFPATH a path containing a number of Amber, Charmm, Cgenff,
Opls, and water forcefields available for use with the -f command-line
option.

Rules file
----------

Any forcefield that is not a forcefield patch must contain a rules file. The
rules file is a JSON dictionary with the following keys and values. Default
values are used for any unspecified keys.

 - **info** (array of strings): Provides citations and any other pertinent
   information. Default: []

 - **vdw_func** (string): The functional form of the VDW interaction.
   Default: ""

 - **vdw_comb_rule** (string): The combining rule for the VDW interaction; 
   legal values depend on the value of **vdw_func**. Default: ""

 - **exclusions** (integer): The exclusion rule; should be 1, 2, 3 or 4.
   Default: Length of es_scale or lj_scale plus 1, or 4 if both es_scale and
   lj_scale are unspecified.

 - **es_scale** (list of floats): Provides a scale factor for each excluded
   nonbonded electrostatic interaction distance. Default: [0,...,0]

 - **lj_scale** (list of floats): Provides a scale factor for each excluded
   nonbonded Lennard-Jones interaction distance. Default: [0,...,0]

 - **plugins** (list of strings): Indicates what functional forms are present
   in the forcefield. Each plugin processes a functional form by reading one
   or more parameter files from the forcefield directory and generating
   interaction terms in the system. Default: []
 
 - **fatal** (bool): Whether ``viparr`` should exit in error if the parameters
   for a required interaction term are not found in the forcefield. Default:
   True
 
 - **nbfix_identifier** (str): A unique identifier label shared by a group of
   forcefields that should be processed together when applying NB-fix
   parameters. Default: ""

Template files
--------------

Any file of name "template..." is treated as a template file. The template file
is a JSON dictionary whose keys are template names and values are dictionaries
describing the atoms, bonds, and other structures present in the template. The
bond topology described by a template must be connected. Atoms in the template
must have unique names; bonds, impropers, etc. may refer to atoms outside the
template by using the names "$1", "$2", etc.

Note that viparr does not employ the concept of "patch" residues; if a residue
has several different chemical forms corresponding to different termini,
phosophorylation states, etc., there should a different template for each form.

Each template may have the following keys and values:

 - **atoms**: ``[ "name", anum, charge, ["btype", "nbtype"]]``.  
   Alternatively: ``[ "name", anum, charge, ["btype"]]```

   *name* is the unique name for the atom.  *anum* is the atomic number.
   *charge* is the partial charge. *btype* is the atom type for bonded
   interactions. If *nbtype* is given, it is used as the atom type for nonbonded
   interactions; otherwise the value of *btype* is used.

 - **bonds**: ``[ "name1", "name2" ]``

   All possible angle and dihedral terms will be autogenerated by ``viparr``
   using the bond topology.

 - **impropers**: ``[ "name1", "name2", "name3", "name4" ]``

 - **cmap**: ``[ "name1", "name2", "name3", "name4", "name5", "name6, "name7", "name8" ]``

 - **pseudos**: ``[ "name", charge, ["btype", "nbtype"], "function", "site_1", ..., "site_n", "pset"``].
   Alternatively, ``[ "name", charge, ["btype"], "function", "site_1", ..., "site_n", "pset"]``.  

   Declares a pseudo-atom for the residue. *charge*, *btype*, and *nbtype* are
   as for **atoms**. *function* is the category of pseudo atom; it may indicate
   the type of virtual site used to position the pseudo atom, and pseudo atoms
   with different functions are processed different by forcefield plugins.
   *site_1* to *site_n* are names of atoms that are used to position the pseudo
   atom. Finally, *pset* is used to distinguish otherwise identical pseudo
   atoms, if the identities of the *sites* are not sufficient.

Parameter files
---------------

``viparr`` expects parameters for a forcefield to be organzized into separate
files based on the functional form of the interaction. Parameter files are JSON
lists of dictionaries; each dictionary corresponds to one parameter set and may
have the following keys:

 - **type** (string or list of strings): Which atom types to match to these
   parameters. The number of types depends on the interaction; there should be
   the same number of types in each parameter set. Types can be delineated as a
   list of strings or ' '-concatenated as a single string.

 - **params** (dictionary): The parameters, with keys being the parameter names
   and values being the corresponding parameter values. Every **params**
   dictionary in a parameter file should have the same keys (parameter names).

Cmap files
----------

For charmm-style cmap terms, parameters may be specified in a forcefield file
called 'cmap'. This is a JSON list of cmap tables, where each cmap table is a
list of [*phi*, *psi*, *energy*] triples specifying the cmap parameters.
The order of cmap tables in this file is important, as they are referred to by
index in the "cmapid" parameter of the 'torsiontorsion_cmap' parameter file.

Extending ``viparr`` functionality
==================================

From C++
--------

The ``viparr`` executable checks the system path VIPARR_PLUGIN_PATH for any .so
shared libraries and dynamically loads them at the start of ``viparr``
execution. This is also performed when the ``viparr`` Python library is imported
into Python. These dynamically-loaded shared libraries should use the DESRES
``plugger`` module to define static initializers that are to be executed when
the libraries are loaded; these initializers may in turn call the static
``Forcefield::RegisterPlugin``, ``Rules::RegisterVDWFunc``, and
``Rules::RegisterVDWCombRule`` functions to allow ``viparr`` to execute with
user-defined plugins and VDW functional forms/combine rules. Please refer to the
``plugger`` library and viparr C++ source code for details.

From Python
-----------

Python modules located in the system path VIPARR_PYPLUGIN_PATH are loaded when
the ``viparr`` Python library is imported into Python. Examples are given below
on using Python to define and add new forcefield plugins or new VDW functional
forms/combine rules. Note that support is not provided to load Python modules
during the execution of the ``viparr`` command-line executable, but ``viparr``
may be executed from within Python, after plugins are loaded, by calling the
``ExecuteViparr`` function in the Python API.

Adding/modifying a forcefield plugin
------------------------------------

Broadly speaking, forcefield plugin can be any function that takes a forcefield
and a system to be parametrized by that forcefield and modifies the system in
some way. The simplest plugin may take a particular list of typified atom
tuples, match them to a particular table of parameters contained in the
forcefield, and add these tuples and matched parameters to the chemical system.
The following example implements a user-defined ``bonds`` plugin, replaces the
default ``bonds`` plugin with this user-defined plugin, and runs ``viparr``.

.. literalinclude:: code/stretch_harm_example.py

Adding/modifying a VDW functional form and combine rule
-------------------------------------------------------

``viparr`` maintains a registry of recognized VDW functional forms and combine
rules; its standard VDW and LJ pairs plugins use the functional form and combine
rule to create the appropriate nonbonded and pairs tables in the system and to
compute the scaled LJ pair interaction terms. The following example defines an
lj12_6 functional form and geometric combine rule, replaces the built-in
versions with these user-defined versions, and runs ``viparr``.

.. literalinclude:: code/lj12_6_example.py
