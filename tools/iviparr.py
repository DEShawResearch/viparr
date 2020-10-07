"""
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
"""

import sys, os, argparse
import msys
import viparr

import subprocess, tempfile, shutil

def configure(parser):
    parser.add_argument("input", help="input dms file")
    parser.add_argument("output", help="output ff dir")
    parser.add_argument("--ffdir", "-d",
                        help="explicit forcefield directory")
    parser.add_argument("--ffname" ,"-f",
                        help="forcefield directory in VIPARR_FFPATH")
    parser.add_argument("--selection" ,"-s", default='all',
                        help="selection of atoms to parametrize (Default: all)")
    parser.add_argument('-t', "--templateonly", action="store_true",
                        help="Do not generate forcefield param tables")
    parser.add_argument("--dontfail", action="store_true",
                        help="Ignore differences between input system and system paramterized with generated forcefield")

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
    configure(parser)
    args = parser.parse_args()
    ff_name = args.ffdir or args.ffname
    if not ff_name:
        parser.error("Provide one of '--ffname' or '--ffdir'")
    ff_path = viparr.find_forcefield(ff_name)
    
    mol = msys.Load(args.input)
    ids = mol.selectIds(args.selection)

    ff = viparr.ImportForcefield(ff_path)
    viparr.Forcefield.ClearParamTables()

    ff_out = viparr._viparr.ExecuteIviparr(mol.asCapsule(), ids, ff._Forcefield.rules(), args.templateonly)

    # store provenance in info section of rules file
    prov = msys._msys.Provenance.fromArgs(sys.argv)
    info = ["Forcefield generated using iviparr",
            "Viparr version: %s" % viparr.version,
            "Timestamp: %s" % prov.timestamp,
            "User: %s" % prov.user,
            "Workdir: %s" % prov.workdir,
            "Cmdline: %s" % prov.cmdline,
            ]
    ff_out.rules().info = info
    viparr._viparr.ExportForcefield(ff_out, args.output);
    if args.templateonly:
        os.unlink(os.path.join(args.output, "rules"))
        return

    # make sure it worked
    with tempfile.NamedTemporaryFile(suffix='.dms') as tmp:
        prog = os.path.dirname(sys.argv[0]) + '/viparr'
        subprocess.check_call([prog, args.input, '-d', args.output, tmp.name])
        rc = subprocess.call(["dms-diff", args.input, tmp.name])
        if rc != 0:
            print("ERROR: dms-diff reports that iviparr did not successfully invert the forcefield:", file=sys.stderr) 
            print("ERROR: A round-trip using viparr did not reproduce the input dms", file=sys.stderr)       
            if not args.dontfail:
                sys.exit(rc)

if __name__=="__main__":
    main()

