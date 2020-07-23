def parse_args():
    from argparse import ArgumentParser
    import os

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("psf", type=os.path.realpath, help="charmm psf file")
    parser.add_argument(
        "ffdir",
        type=os.path.realpath,
        help="Output template file into a viparr forcefied",
    )
    parser.add_argument(
        "-f",
        "--fuse",
        default=False,
        action="store_true",
        help="Fuse all residues into a single template",
    )
    parser.add_argument(
        "-s",
        "--selection",
        default="all",
        type=str,
        help="msys atom selection of atoms in the psf to convert into templates",
    )
    return parser.parse_args()


def psftemplate(psfmol, fuse, selection):
    import viparr
    from viparr import ffconverter

    viparrff = ffconverter.initializeSimpleForcefield(
        [0, 0, 1.0], [0, 0, 1.0], "charmm"
    )

    mysel = selection
    if mysel != "all" and "same" not in mysel:
        mysel = "same residue as "+mysel
    only = psfmol.select(mysel)

    nbtable = psfmol.table("nonbonded")

    reslist = psfmol.updateFragids() if fuse else psfmol.residues
    residues = []
    for r in reslist:
        atoms = [a for a in r if a in only]
        if len(atoms) == 0 : continue
        residues.append(atoms)

    if fuse:
        for frag in residues:
           frag[0].residue.name = "_".join(sorted({a.residue.name for a in frag}))

    nbtable = psfmol.table("nonbonded")
    for ires, resatoms in enumerate(residues):
        tsys = viparr.TemplatedSystem()
        res = tsys.system.addResidue()
        res.name = resatoms[0].residue.name

        qsum = 0.0
        amap = {}
        anames = {}
        for atom in resatoms:
            a = res.addAtom()
            a.name = atom.name
            assert a.name not in anames, "duplicated anames detected"
            a.atomic_number = atom.atomic_number
            assert a.atomic_number > 0, "Virtuals not supported"
            a.charge = atom.charge
            qsum += a.charge
            atypes = nbtable.findWithAll([atom])
            assert (
                len(atypes) == 1
            ), f"There should only be one atype found. Got {atypes}"
            atype = atypes[0].param["type"]
            tsys.setTypes(a, atype, atype)
            amap[atom] = a

        dangling = {
            other
            for atom in resatoms
            for other in atom.bonded_atoms
            if other not in resatoms
        }
        for i, atom in enumerate(dangling):
            a = res.addAtom()
            a.name = f"${i}"
            a.atomic_number = -1
            amap[atom] = a

        resbonds = {b for a in resatoms for b in a.bonds}
        for b in resbonds:
            a0 = amap[b.first]
            a1 = amap[b.second]
            a0.addBond(a1)

        if "improper_harm" in psfmol.table_names:
            for impr in psfmol.table("improper_harm").terms:
                alist = [amap[a] for a in impr.atoms if a in amap]
                asum = sum([1 if a.atomic_number > 0 else 0 for a in alist])
                if len(alist) != 4 or asum < 3:
                    continue
                tsys.addImproper(alist)

        if "torsiontorsion_cmap" in psfmol.table_names:
            for cmap in psfmol.table("torsiontorsion_cmap").terms:
                alist = [amap[a] for a in cmap.atoms if a in amap]
                asum = sum([1 if a.atomic_number > 0 else 0 for a in alist])
                if len(alist) != 8 or asum < 6:
                    continue
                print(alist)
                tsys.addCmap(alist)

        print(res.name, qsum)
        ffconverter.addTemplateData(viparrff.typer, tsys, True)
    return viparrff


def main():
    import msys
    import viparr

    args = parse_args()
    mol = msys.Load(args.psf)

    ff = psftemplate(mol, args.fuse, args.selection)
    viparr.ExportForcefield(ff, args.ffdir)
