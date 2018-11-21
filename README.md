Viparr is a set of tools and libraries for building forcefields
for molecular systems suitable for molecular dynamics simulations.
It requires the msys library, distributed by D.E. Shaw Research.
Forcefield files for commonly used forcefields may be found in the
D. E. Shaw Research `viparr-ffpublic` github repository.

Building Viparr
---------------

Building viparr requires that you have first built and installed msys; see the associated instructions.  The lpsolve and InChI optional features are not required for viparr.  Like msys, viparr requires Scons, a C++ build tool; numpy headers for each python version you want to support; and boost headers and libraries (`filesystem`, `system`, and `program_options`).  

To build and install viparr, use the scons tool:

    PYTHONPATH=external:$PYTHONPATH scons -j4 PYTHONVER=37 PREFIX=$HOME/local install

where `PREFIX` is the location of the msys installation.

Documentation for viparr is in sphinx format.  To build the documentation, first install viparr, then run the following in the `doc` directory:

    PYTHONPATH=$HOME/local/lib/python make clean html


Running Viparr
--------------

Fetch the forcefield files from the `viparr-ffpublic` repository, and set the environment variable `VIPARR_FFPATH` to the full path of the `ff` directory located within the `viparr-ffpublic` checkout.  To parameterize, for example, a molecular system containing protein, lipids, water, and ions, using a set of Amber forcefields, specify the name of each forcefield directory as an argument to `viparr`:

    viparr input.dms output.dms -f aa.amber.ff99SB-ILDN -f lipid.amber.lipid17 -f water.tip3p

If you have a viparr-compatible forcefield directory that isn't contained within `viparr-ffpublic`, you may specify it directly with the `-d` option.

Viparr is distributed with tools to convert forcefields from CHARMM and Amber formats.


