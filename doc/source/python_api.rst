``````````
Python API
``````````

The ``viparr`` Python library provides Python wrappers around the core classes
and functionality of the ``viparr`` executable (written in C++). Its purpose is
to provide forcefield developers with in-memory representations of Forcefield
objects and modularized access to the individual subroutines (atom typing,
parameter matching, etc.) of ``viparr`` execution.

Examples
========

Constructing a new template
------------------------------------

As an example, the Python API may be used to help in the constructon of a new 
forcefield template for a given molecule:

.. literalinclude:: code/construct_template_example.py


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


API Reference
=============

.. automodule:: viparr
    :members:
