````````````````````````
Other command-line tools
````````````````````````

build_constraints
=================

::

  $ build_constraints system.dms [-s selection] [-k] [-x type]*

Builds or rebuilds constraints for hydrogen groups in a chemical system. This is
performed by default at the end of ``viparr`` execution. Each hydrogen group
consists of a heavy atom and all hydrogens bonded to it. A heavy atom with `n`
hydrogens will be constrained with a constraint of type `ahn`. In the case of
water, the constraint type will be `hoh`.

The constraint length parameters (and angle parameter, in the case of water)
are obtained from the ``stretch_harm`` and ``angle_harm`` tables. Stretch and
angle terms that are constrained by a constraint group are marked by setting the
`constrained` field in the ``stretch_harm`` or ``angle_harm`` table to 1, unless
the `-k` option is given. Constrained terms can safely be ignored by simulation
programs as long as the constraints are evaluated.

Options:

.. program:: build_constraints

.. cmdoption:: -s selection

   Selection of atoms; must correspond to a subset of complete connected molecules. Default is 'all'.

.. cmdoption:: -k

   Do not set the 'constrained' field in stretch-harm and angle-harm tables to true.

.. cmdoption:: -x type

   Do not build constraints of a given type; type can be 'hoh', 'ah1', 'ah2', etc.

construct_template
==================

::

  $ construct_template input_system selection output_template [--by-residue]

Creates a forcefield template from a selection of atoms in a chemical system.
This populates only the atom and bond lists of the template using the bond
topology of the system; it does not create lists of impropers, exclusions, etc.
To create a full template from a parametrized chemical system, use ``iviparr``.

Options:

.. program:: construct_template

.. cmdoption:: --by-residue

   Whether to create a separate template for each residue in the ``selection``.

dms_from_template
=================

::

  $ dms_from_template tplfile [tplname]*

.. program:: dms_from_template

Constructs DMS system files from templates. This will construct one DMS file
for each specified template and export it to `tplname.dms`. The DMS file will
not contain forcefield information.

draw_template
=============

Generates a PNG image of a forcefield template.

::

  $ draw_template template_file template_name [--dotfile dotfile_path] \
      [--pngfile pngfile_path] [--nodisplay]

Creates a PNG image of a forcefield template, specified by a template file and
the name of a particular template in that file. The program constructs a
dot-format text file of the template graph, saves it to the given output path,
calls the shell command `neato -Tpng dotfile -o pngfile` to convert the dot file
to a PNG image file using graphviz, saves the image file to the given output
path, and displays the PNG image using gthumb.

Options:

.. program:: draw_template

.. cmdoption:: --dotfile dotfile_path

   Location to save output dot file. Default: `topology.dot`

.. cmdoption:: --pngfile pngfile_path

   Location to save output png file. Default: `topology.png`

.. cmdoption:: --nodisplay

   Do not display graphviz output in gthumb

iviparr
=======

::

  $ iviparr input_system output_ff [-s selection] \
      [-f ffname | -d ffdir] [--templateonly]

Creates a forcefield from a parametrized chemical system. The selected atoms
should correspond to a set of fragments that were parametrized using a single
forcefield. A separate template is created for each residue of the selected
atoms, and parameter tables are created from the forcefield parameters in the
system. As certain aspects of the rules file cannot be deduced from the chemical
system, the program requires an input forcefield, specified by -f or -d, to use
only the rules file contained in that forcefield. ``iviparr`` will exit in error
if it is unable to construct a template-based forcefield that can yield the
given parametrized system when run using ``viparr``, e.g. if there are multiple
residues of the system with identical bond topology but different parameters.

Options:

.. program:: iviparr

.. cmdoption:: -s selection

   Selection of atoms for which to create forcefield. Default: all.

.. cmdoption:: -f ffname, -d ffdir

   A forcefield with a rules file that matches the one used to parametrize the selected atoms.

.. cmdoption:: --templateonly

   Create forcefield with rules and templates only (no parameter files).

merge_forcefields
=================

::

  $ merge_forcefields [(-f ffname | -d ffdir) [-m ffpatch | -a ffpatch]* ] output_ff

Merges forcefield patches into a source forcefield. Using the ``output_ff`` to
parametrize a system with ``viparr`` is equivalent to directly specifying the
merge in the ``viparr`` call. Forcefield patches are merged sequentially into
the source forcefield in the order in which they are specified.

Plugins are merged in the rules files with the plugins of the patch forcefield
coming after those of the source forcefield. VDW functional forms, VDW combine
rules, exclusions rules, and scale factors, if specified in both the source and
patch forcefields, must agree. Templates in the patch forcefield are added to 
the source forcefield and overwrite those in the source forcefield of the same
name. Parameters in the patch forcefield overwrite those in the source
forcefield having the same `type` fields (i.e. matching the same atom types);
the merged parameter file will begin with the parameters in the patch forcefield
not found in the source forcefield, followed by the parameters of the source
forcefield possibly overwritten by ones from the patch matching the same types.
Cmap tables in the patch forcefield, if present, overwrite all cmap tables in
the source forcefield.

Options:

.. program:: merge_forcefields

.. cmdoption:: -f ffname, -d ffdir

   The source forcefield.

.. cmdoption:: -m ffpatch

   A forcefield patch to be merged into the source forcefield.

.. cmdoption:: -a ffpatch

   Like -m ffpatch, except that elements of ffpatch cannot overwrite those of the source forcefield.

optimize_vsitedefs
==================

::

  $ optimize_vsitedefs system.dms

.. program:: optimize_vsitedefs

This program tries to turn lcXn virtual site definitions into lcX definitions
by checking if the defining atoms for the virtual site depend completely on
constrained bonds. This routine is run automatically whenever constraints are
built either by ``viparr`` or ``build_constraints``.

viparr_neutralize
=================

::

  $ viparr_neutralize input_system output.dms [-p cation] [-n anion] [-c chain] \
      [-C chain2] [-s solute_pad] [-i ion_pad] [-m concentration] \
      [[--ffname ff_name | --ffdir ff_dir] | --no-ff] [--random-seed seed]

Adds or removes ions from a chemical system in order to neutralize the charge
and attain a specified ionic concentration. A forcefield may be provided to
parametrize any added ions. Randomly selected water molecules may be removed
from the system to make room for the added ions. In the following, counterions
are those with opposite charge from the solute (required to neutralize the
charge of the system), and counter-counterions are those of the same charge as
the solute (possibly required to achieve a desired ionic concentration).

Options:

.. program:: viparr_neutralize

.. cmdoption:: -p cation

   The species of cation to add or remove; may be NA or K. Default: NA.

.. cmdoption:: -n anion

   The species of anion to add or remove; may be CL. Default: CL.

.. cmdoption:: -c chain

   The chain name for added counterions.

.. cmdoption:: -C chain2

   The chain name for added counter-counterions.

.. cmdoption:: -s solute_pad

   Minimum distance between added ions and any non-water molecule in the system.

.. cmdoption:: -i ion_pad

   Minimum distance between two added ions.

.. cmdoption:: -m concentration

   Molar concentration of counter-counterions. The number of such ions will be
   given by int((concentration / 55.345) * (nwater - nions)) where nwater is
   the total number of waters in the original system and nions is the number of
   counterions needed to achieve charge neutrality.

.. cmdoption:: --ffname ff_name, --ffdir ff_dir

   A forcefield to be applied to any added ions, specified by name or absolute path.

.. cmdoption:: --no-ff

   No forcefield is to be applied; existing forcefield in the system is removed.
   A required option if no forcefield option is provided.

.. cmdoption:: --random-seed seed

   Seed for the random number generator to determine which water molecules are replaced.

viparr_solvate
==============

::

  $ viparr_solvate sys1.dms [sys2.dms] [sys3.dms] [-d dims] [-c center] \
      [-n chain] [--ffname ff_name | --ffdir ff_dir] \
      [--without-constraints] [--verbose]

Adds water or molecules of a user-specified solvent to a system. If one DMS
argument is supplied, it is interpreted as `watbox_out.dms` and a system
containing only water is saved to this file. If two DMS arguments are supplied,
they are interpreted as `solute_in.dms` and `solution_out.dms`; water molecules
are tiled around the system in `solute_in.dms` and the resulting system is saved
to `solution_out.dms`. If three DMS arguments are supplied, they are interpreted
as `solute_in.dms`, `watbox_in.dms`, and `solution_out.dms`; the solvent box
given in `watbox_in.dms` is tiled around `solute_in.dms` and the resulting
system is saved to `solution_out.dms`.

A forcefield may be specified to parametrize the default water box or the
user-provided `watbox_in.dms` system; the forcefield will not be applied to the
solute system. If a forcefield is applied, then constraints are built by default
at the end of execution. If no forcefield is applied, then the default water box
is tiled unparametrized in the case of one or two DMS arguments, while a
user-provided `watbox_in.dms` solvent keeps its existing forcefield parameters
in the case of three DMS arguments.

.. program:: viparr_solvate

.. cmdoption:: -d dims

   Dimensions of the output system, specified as 3 comma-separated values or a
   single float value for a square box. Default takes the maximum of the
   dimensions of the input solute system and the dimensions of the input solvent
   box.

.. cmdoption:: -c center

   Center of the output system as 3 comma-separated values. Default: 0,0,0.

.. cmdoption:: -n chain

   The chain name for added water or solvent molecules.

.. cmdoption:: --ffname ff_name, --ffdir ff_dir

   A forcefield to be applied to any added water or solvent molecules.

.. cmdoption:: --without-constraints

   Do not build constraints for the added water or solvent molecules.
