def parse_args():
    from argparse import ArgumentParser
    import os

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("psf", type=os.path.realpath, help="charmm psf file")
    parser.add_argument(
        "ffdir", type=os.path.realpath, help="Output template file into a viparr forcefied"
    )
    parser.add_argument(
        "-f",
        "--fuse",
        default=False,
        action="store_true",
        help="Fuse all residues into a single template",
    )
    return parser.parse_args()


def psftemplate(psfmol, fuse):
    import viparr
    from viparr import ffconverter

    viparrff = ffconverter.initializeSimpleForcefield(
        [0, 0, 1.0], [0, 0, 1.0], "charmm"
    )

    nbtable = psfmol.table("nonbonded")

    if (fuse):
        residues = psfmol.updateFragids()
        for frag in residues:
            name = []
            for res in sorted({(a.residue.id, a.residue) for a in frag}):
                name.append(res[1].name)
            frag[0].residue.name = "_".join(name)
    else:
        residues = [r.atoms for r in psfmol.residues]

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
            a.name = f"$i"
            a.atomic_number = -1
            amap[atom] = a

        resbonds = {b for a in resatoms for b in a.bonds}
        for b in resbonds:
            a0 = amap[b.first]
            a1 = amap[b.second]
            a0.addBond(a1)

        if "improper_harm" in psfmol.table_names:
            for impr in psfmol.table("improper_harm").terms:
                print([a.fullname for a in impr.atoms])
                alist = [amap[a] for a in impr.atoms if a in amap]
                if alist != 4 or sum([1 if a.atomic_number > 0 else 0 for a in alist]) <3:
                    continue
                tsys.addImproper(alist)

        if "torsion_torsion_cmap" in psfmol.table_names:
            for cmap in psfmol.table("torsion_torsion_cmap").terms:
                alist = [amap[a] for a in cmap.atoms if a in amap]
                if alist != 8 or sum([1 if a.atomic_number > 0 else 0 for a in alist]) <7:
                    continue
                tsys.addCmap(alist)

        print(res.name, qsum)
        ffconverter.addTemplateData(viparrff.typer, tsys, True)
    return viparrff


def main():
    import msys
    import viparr

    args = parse_args()
    mol = msys.Load(args.psf)

    ff = psftemplate(mol, args.fuse)
    viparr.ExportForcefield(ff, args.ffdir)