
"""
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
"""

import msys
import viparr

import argparse
import os
import shutil
import subprocess
import tempfile


class PrintAvailableForcefields:
    def __mod__(self, prog):
        usage_lines = ["viparr input.dms output.dms [ options ]"]
        usage_lines.append("VIPARR version {}".format(viparr.__version__))
        ffpath = viparr.get_ffpath()
        if not ffpath:
            usage_lines.append("No available forcefields; set VIPARR_FFPATH " + \
                "by loading module viparr-ff/*/data")
        else:
            ffs = viparr.list_forcefields(ffpath)
            usage_lines += ["VIPARR_FFPATH: %s" % ffpath]
            usage_lines += ["All available forcefields", "--------------------------"]
            fmt = "%-"+str(max(map(len,ffs)))+"s --- %s"
            for name, path in sorted(ffs.items()):
                rules_file = '%s/rules' % path
                rules = viparr.ImportRules(rules_file)
                info = rules.info[0] if rules.info else "No information"
                usage_lines.append(fmt % (name, info))
        return "\n".join(usage_lines)

class FFDirAction(argparse.Action):
    def __call__(self, parser, namespace, value, option_string):
        if not os.path.isdir(value):
            raise argparse.ArgumentError(self, "'%s' not found or is not a directory" % value)
        print('Importing forcefield from %s' % value)
        ff = viparr.ImportForcefield(value)
        namespace.fflist.append(ff)

class FFNameAction(argparse.Action):
    def __call__(self, parser, namespace, value, option_string):
        src = viparr.find_forcefield(value)
        print('Importing forcefield from %s' % src)
        namespace.fflist.append(viparr.ImportForcefield(src))

class FFMergeAction(argparse.Action):
    def __call__(self, parser, namespace, value, option_string):
        if not namespace.fflist:
            raise argparse.ArgumentError(self, "No previous forcefield specified with -d/-f")
        target = namespace.fflist[-1]

        src = viparr.find_forcefield(value)
        print("Importing forcefield patch from %s" % src)
        patch = viparr.ImportForcefield(src, False);
        viparr.MergeForcefields(target, patch)

class FFAppendAction(argparse.Action):
    def __call__(self, parser, namespace, value, option_string):
        if not namespace.fflist:
            raise argparse.ArgumentError(self, "No previous forcefield specified with -d/-f")
        target = namespace.fflist[-1]
        print("Importing append-only forcefield patch from %s" % value)
        patch = viparr.ImportForcefield(value, False);
        viparr.MergeForcefields(target, patch, True)

def parser():
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            usage=PrintAvailableForcefields())
    parser.add_argument("input", help="input dms file")
    parser.add_argument("output", help="output dms file")
    parser.add_argument("--ffdir", "-d", dest='fflist', default=[], action=FFDirAction,
                        help="explicit forcefield directory")
    parser.add_argument("--ffname" ,"-f", dest='fflist', action=FFNameAction,
                        help="forcefield directory in VIPARR_FFPATH")
    parser.add_argument("--merge", "-m", dest='fflist', action=FFMergeAction,
                        help="merge forcefield patch; may overwrite existing forcefield components")
    parser.add_argument("--append", "-a", dest='fflist', action=FFAppendAction,
                        help="merge forcefield patch; only append to existing forcefield components")
    parser.add_argument("--ffde",
                        help="path to FFDE forcefield")
    parser.add_argument('--ligand-files', nargs='*',
                        help="one or more parameterized ligand dms files")
    parser.add_argument('--ligand-selection', default='none',
                        help="selection for ligand atoms")
    parser.add_argument('--match-ligand-bond-stereo', action='store_true',
                        help="Require bond stereo match in ligand input files")
    parser.add_argument('--match-ligand-tet-stereo', action='store_true',
                        help="Require tetrahedral stereo match in ligand input files")
    parser.add_argument('--match-ligand-hydrogen', action='store_true',
                        help="Require inchi hydrogen layer match in ligand input files")
    parser.add_argument('--exhaustive-ligand-matching', action='store_true',
                        help="Search all possible mappings of ligand inputs")
    parser.add_argument("--selection" ,"-s", default='all',
                        help="selection of atoms to parametrize (Default: all)")
    parser.add_argument("--reorder-ids", action="store_true", help="reorder atom ids so that pseudos are next to parents")
    parser.add_argument("--rename-atoms", action="store_true", help="copy template atom names to system")
    parser.add_argument("--rename-residues", action="store_true", help="copy template residue names to system")
    parser.add_argument("--without-constraints", action="store_false", dest='with_constraints', help="do not build constraints")
    parser.add_argument("--without-fix-masses", action="store_false", dest='fix_masses',
                        help="do not equate masses for atoms of the same element")
    parser.add_argument("--verbose-plugins", action="store_true", help="print debug messages from plugin load")
    parser.add_argument("--non-fatal", action="store_true", help="do not exit in error if a required term is unmatched")
    parser.add_argument("--verbose-matching", action="store_true", help="print the template matched by each residue")

    parser.add_argument("--make-rigid", action="store_true",
                                help="Replaces constraint_ah1, _ah2, and _ah3 constraints with alternative " \
                                "constraint_ah1R, _ah2R, and _ah3R constraints.")
    return parser

def run_viparr(args):
    structure_only = True if args.selection == 'all' else False
    if args.ligand_files and args.ligand_selection=='none':
        raise RuntimeError("Missing --ligand-selection")
    print("Importing structure from %s" % args.input)
    if args.ffde:
        if args.fflist:
            raise RuntimeError("Cannot specify both ffpath and other forcefield options")
        print("Exporting parametrized system to %s" % args.output)
        viparr.ExecuteFFDETypify(args.input, args.output, args.ffde)
        return

    mol = msys.Load(args.input, structure_only=structure_only)
    ffs = [ff._Forcefield for ff in args.fflist]
    ids = mol.selectIds('(%s) and not (%s)' % (args.selection, args.ligand_selection))
    if not ids:
        if args.ligand_selection == 'none':
            raise RuntimeError("no atoms selected and no ligands to parameterize")
        print("No atoms to parameterize using forcefields")
    else:
        compile_plugins = True  # ??
        optimize_vsite_defs = True
        viparr._viparr.ExecuteViparr(mol.asCapsule(), ffs, ids,
                args.rename_atoms, args.rename_residues, args.with_constraints,
                optimize_vsite_defs,
                args.fix_masses, not args.non_fatal,
                compile_plugins, args.verbose_plugins, args.verbose_matching)

    if args.ligand_files:
        ligands = [msys.Load(f) for f in args.ligand_files]
        mol = viparr.ApplyLigandForcefields(mol, ligands, args.ligand_selection,
                args.rename_atoms, args.rename_residues,
                verbose=args.verbose_matching,
                match_bond_stereo=args.match_ligand_bond_stereo,
                match_tet_stereo=args.match_ligand_tet_stereo,
                match_hydrogen=args.match_ligand_hydrogen,
                exhaustive_matching=args.exhaustive_ligand_matching)



    if args.reorder_ids:
        print("Reordering pseudo IDs")
        mol = viparr.ReorderIDs(mol)

    mol.coalesceTables()
    mol = mol.clone()

    # FIXME: do the post-processing steps before saving
    print("Exporting parametrized system to %s" % args.output)
    msys.Save(mol, args.output)

def main():
    args = parser().parse_args()
    run_viparr(args)

    # FIXME: get this done without an extra read/write step.
    if args.make_rigid:
        actual_out = os.path.abspath(args.output)
        tmp_out = tempfile.mktemp()
        shutil.move(actual_out, tmp_out)
        subprocess.check_call(["viparr-make-rigid", tmp_out, actual_out])

    print("VIPARR exited successfully")

# vim: filetype=python
