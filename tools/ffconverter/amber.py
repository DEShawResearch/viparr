from __future__ import print_function

import os,numpy
import msys, viparr
from . import helper as viparrConversionHelper

# Do not change these unless you know what you are doing.
# You will be in for a world of hurt if you change these and then attempt
# to merge or append forcefields together in viparr
amberEsScales=[0,0,0.8333]
amberVdwScales=[0.0, 0.0, 0.5]

def convertPrmtopToDMS(prmtop,crd):
    if(not (os.path.exists(prmtop) and os.path.exists(crd))):
        raise RuntimeError("Amber prmtop or crd file missing")
    print('Converting AMBER files to dms')
    molNew = msys.LoadPrmTop(prmtop, False)
    msys.ReadCrdCoordinates(molNew, crd)
    return molNew

def createViparrForcefieldFromAmberDMS(molFF, symQ=True, ff=None):
    from collections import defaultdict
    if(ff is None):
        viparrff=viparrConversionHelper.initializeSimpleForcefield(amberEsScales, amberVdwScales)
    else:
        viparrff=ff

    # extract the atomtype data first so we can fill in the parameter files and make a template with atomtypes
    atypes=[]
    nbtable = molFF.getTable('nonbonded')
    for i, t in enumerate(nbtable.terms):
        if(len(t.atoms) != 1 or t.atoms[0].id != i):
            raise RuntimeError("Consistency Check Failed: %d!=1 or %d!=%d" %(len(t.atoms), t.atoms[0].id, i))
        pdict={k : t.param[k] for k in list(t.param.keys())}
        viparrConversionHelper.addOrUpdateParameterData(viparrff, "vdw1", pdict, False, False)
        atypes.append(t.param["type"].replace('*', '&'))
    if(len(atypes) != molFF.natoms):
        raise RuntimeError("Consistency Check Failed. Missing atomstypes: %d != %d"%(len(atypes), molFF.natoms))

    acharges=[]
    tsys = viparr.TemplatedSystem()
    res = tsys.system.addResidue()
    res.name = molFF.name
    for atom, atype in zip(molFF.atoms, atypes):
        if(atom.atomic_number < 1):
            raise RuntimeError("non-real atoms are not supported")
        a = res.addAtom()
        a.name = atom.name
        a.atomic_number = atom.atomic_number
        a.mass = atom.mass
        a.charge = atom.charge
        tsys.setTypes(a, atype, atype)
        acharges.append(atom.charge)
        pdict={"type" : atype, "amu" : atom.mass, "memo" : ""}
        viparrConversionHelper.addOrUpdateParameterData(viparrff, "mass", pdict, False, False)
    for bond in molFF.bonds:
        res.system.atom(bond.first.id).addBond(res.system.atom(bond.second.id))
    viparrConversionHelper.addTemplateData(viparrff.typer, tsys,symQ)

    bads = viparr.GetBondsAnglesDihedrals(molFF)
    # This shouldnt happen it there were no fake atoms above
    assert(len(bads["pseudo_bonds"]) == 0)
    non_pseudo_bonds = [ (molFF.atom(i), molFF.atom(j)) for i, j in bads["non_pseudo_bonds"] ]
    angles = [ (molFF.atom(i), molFF.atom(j), molFF.atom(k)) for i, j, k in bads["angles"] ]
    dihedrals = [ (molFF.atom(i), molFF.atom(j), molFF.atom(k), molFF.atom(l)) for i, j, k, l in bads["dihedrals"] ]

    termCounts=defaultdict(int)
    for table in molFF.tables:
        if(table.name == "nonbonded"): continue
        if(table.name in ["stretch_harm", "angle_harm", "dihedral_trig",
                          "improper_trig", "improper_harm", "improper_anharm"]):
            plist=[]
            for t in table.terms:
                ttype = [atypes[a.id] for a in t.atoms]
                pdict = {k : t.param[k] for k in list(t.param.keys())}
                atoms = t.atoms
                if(not table.name.startswith("improper") and atoms[-1].id < atoms[0].id):
                    atoms = atoms[::-1]
                if(table.name == 'dihedral_trig'):
                    if(tuple(atoms) not in dihedrals):
                        # This is really an improper, use original atom order
                        plist.append(("improper_trig", t.atoms, pdict))
                        continue
                    if(len(plist) and plist[-1][1] == atoms):
                        plist[-1][2].append(pdict)
                    else:
                        plist.append((table.name, atoms, [pdict]))
                else:
                    plist.append((table.name, atoms, pdict))

            for tname, atoms, pdicts in plist:
                ttype=[atypes[a.id] for a in atoms]
                if(tname.startswith("improper")):
                    tsys.addImproper([tsys.system.atom(a.id) for a in atoms])
                elif(ttype[-1] < ttype[0]):
                    ttype = ttype[::-1]
                ttype=" ".join(ttype)

                if(isinstance(pdicts, list) and len(pdicts) == 1):
                    assert(tname == "dihedral_trig")
                    pdicts = pdicts[0]

                if(isinstance(pdicts, list) ):
                    assert(tname == "dihedral_trig")
                    for d in pdicts: d["type"] = ttype
                    viparrConversionHelper.addOrUpdateParameterData(viparrff, tname, pdicts, False, True)
                else:
                    pdicts["type"] = ttype
                    viparrConversionHelper.addOrUpdateParameterData(viparrff, tname, pdicts, False, False)

                termCounts[tname] += 1

        elif(table.name == "exclusion"):
            # exclusions should be the end atoms of the bonds, angles and torsions
            exclRef = set()
            for data in [non_pseudo_bonds, angles, dihedrals]:
                for term in data:
                    ai, aj =term[0], term[-1]
                    if(aj.id < ai.id):
                        exclRef.add((aj, ai))
                    else:
                        exclRef.add((ai, aj))
            for t in table.terms:
                ai, aj = t.atoms
                if(aj.id < ai.id):
                    assert((aj, ai) in exclRef)
                else:
                    assert((ai, aj) in exclRef)
        elif(table.name == "pair_12_6_es"):
            # make sure pair terms utilize the correct (simple) combining and scale factors.
            # would could probably support *slightly* more complicated scale combinations if necessary
            # (like if we convert glycam with full 1-4 interactions)
            # This is a check only, we dont actually do anything with these parameters
            scaleEsExact = 5./6.
            scaleVDWExact = 0.5
            for t in table.terms:
                if(len(t.atoms)!=2):
                    raise RuntimeError("Pair Terms are malformed: %d != 2"%(len(t.atoms)))
                qSigEpsRef = []
                for a in t.atoms:
                    termRef = nbtable.term(a.id)
                    qSigEpsRef.append((acharges[a.id], termRef["sigma"], termRef["epsilon"]))
                qij, aij, bij = t.param["qij"], t.param["aij"], t.param["bij"]

                qijRef = qSigEpsRef[0][0]*qSigEpsRef[1][0]
                qComb0 = scaleEsExact*qijRef
                qComb1 = amberEsScales[-1]*qijRef
                if(not numpy.isclose(qij, qComb0) and not numpy.isclose(qij, qComb1)):
                    delta = qij - qComb0
                    scale = qij / (qComb0 / scaleEsExact)
                    raise RuntimeError("Pair Term qij != esScale*qi*qj: %f != %f (delta=%f apparentScale=%f)"%(qij, qComb0, delta, scale))
                sigComb = pow(0.5 * (qSigEpsRef[0][1] + qSigEpsRef[1][1]), 6)
                epsComb = 4.0 * scaleVDWExact * (qSigEpsRef[0][2] * qSigEpsRef[1][2])**0.5
                bComb = epsComb * sigComb
                aComb = bComb * sigComb
                if(not numpy.isclose(aij, aComb)):
                    delta = aij - aComb
                    scale = aij / (aComb / scaleVDWExact)
                    raise RuntimeError("Pair Term aij != vdwScale*eps_ij*sigma_ij^12: %f != %f (delta=%f apparentScale=%f)"%(aij, aComb, delta, scale))
                if(not numpy.isclose(bij, bComb)):
                    delta = bij - bComb
                    scale = bij / (bComb / scaleVDWExact)
                    raise RuntimeError("Pair Term bij != vdwScale*eps_ij*sigma_ij^6: %f != %f (delta=%f apparentScale=%f)"%(bij, bComb, delta, scale))
        elif(table.name.startswith("constraint_")):
            pass
        else:
            raise RuntimeError("Unsupported table for conversion: %s"%(table.name))
    if( (len(non_pseudo_bonds) != termCounts["stretch_harm"]) or
        (len(angles) != termCounts["angle_harm"]) or
        (len(dihedrals) != termCounts["dihedral_trig"]) ):
        raise RuntimeError("Consistency check failed. There are missing valence terms")

    return viparrff
