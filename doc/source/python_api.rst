``````````
Python API
``````````

The ``viparr`` Python library provides Python wrappers around the core classes
and functionality of the ``viparr`` executable (written in C++). Its purpose is
to provide forcefield developers with in-memory representations of Forcefield
objects and modularized access to the individual subroutines (atom typing,
parameter matching, etc.) of ``viparr`` execution.

Constructing a new template
------------------------------------

As an example, the Python API may be used to help in the constructon of a new 
forcefield template for a given molecule:

.. literalinclude:: code/construct_template_example.py



API Reference
=============

.. automodule:: viparr
    :members:
