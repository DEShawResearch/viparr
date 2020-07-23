``````
Viparr
``````

Viparr is a tool that applies forcefield parameters to a chemical system
file, producing an output file that can be used to perform molecular
dynamics simulations with Anton or Desmond.

It contains a collection of smaller command-line tools to handle a number
of forcefield-related tasks, as well as a Python library enabling easy
customization and extension of certain aspects of viparr's execution.

Quickstart
==========

Before running any viparr tool, your chemical system must be fully
modeled; i.e., all required hydrogens are present, and extraneous
chemical groups (such as those commonly found as co-crystallization
agents) removed.

Viparr works by matching **templates** in a forcefield library to 
**residues** in your chemical system, using just two attributes:
the atomic number of the atoms, and the bond topology.  If no match
is found, viparr will fail to parameterize the system and try to
tell you where the problem is.

Here is a typical viparr invocation to parameterize a system containing
protein, water, and possibly ions::

    garden load viparr-ff/2.1.4c7/data
    viparr system.dms system.ff.dms -f aa.amber.ff99SB-ILDN -f water.tip3p

The forcefields **aa.amber.ff99SB-ILDN** and **water.tip3p**  are supplied
by the ``viparr-ff`` garden module.


To add a water box to a chemical system, use ``viparr_solvate``::

    viparr_solvate system.ff.dms system.solv.dms --ffname water.tip3p -t 10 -c protein

This command tiles a pre-equilibrated water box around the atoms in the
input system, such that there is at least 10A from any system atom to the
edge of the periodic cell, as well as reasonable distances between water
atoms and system atoms.

To neutralize a solvated chemical system, use ``viparr_neutralize``::

    viparr_neutralize system.solv.dms system.neut.dms -m 0.05 --ffname aa.amber.ff99SB-ILDN

This command replaces water molecules in the system with positive and
negative ions in order to reach the specified concentration (50 mM in
thise case) and to satisify charge neutrality.


.. toctree::
   :maxdepth: 2

   tools
   forcefields
   python_api
   release_notes.rst

