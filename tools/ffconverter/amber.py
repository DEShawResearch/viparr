from __future__ import print_function

import os,numpy
import msys, viparr
from . import helper as viparrConversionHelper

# Do not change these unless you know what you are doing.
# You will be in for a world of hurt if you change these and then attempt
# to merge or append forcefields together in viparr
amberEsScales=[0, 0, 0.8333]
amberVdwScales=[0.0, 0.0, 0.5]

def convertPrmtopToDMS(prmtop,crd):
    if(not (os.path.exists(prmtop) and os.path.exists(crd))):
        raise RuntimeError("Amber prmtop or crd file missing")
    print('Converting AMBER files to dms')
    molNew = msys.LoadPrmTop(prmtop, False)
    msys.ReadCrdCoordinates(molNew, crd)
    return molNew

def createViparrForcefieldFromAmberDMS(molFF, selection, es_scale, vdw_scale, tname=None, symQ=True, ff=None):
    from collections import defaultdict

    my_es_scales = list(amberEsScales)
    my_vdw_scales = list(amberVdwScales)

    if es_scale != my_es_scales[-1]:
        print(f"Warning: You are using a non-standard es scale factor for amber: {es_scale}")
        my_es_scales[-1] = es_scale

    if vdw_scale != my_vdw_scales[-1]:
        print(f"Warning: You are using a non-standard vdw scale factor for amber: {vdw_scale}")
        my_vdw_scales[-1] = vdw_scale

    if(ff is None):
        viparrff = viparrConversionHelper.initializeSimpleForcefield(my_es_scales, my_vdw_scales)
    else:
        viparrff = ff


    selected_ids = set(molFF.selectIds(selection))

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

    external_ids = set([ o.id for a in molFF.atoms for o in a.bonded_atoms if (a.id in selected_ids and o.id not in selected_ids)])
    keep_ids = external_ids | selected_ids

    iexternal = 0
    keep_ids = external_ids | selected_ids
    acharges=[]
    tsys = viparr.TemplatedSystem()
    res = tsys.system.addResidue()
    if tname is None:
        res.name = molFF.name
    else:
        res.name = tname

    idmap = dict()
    for atom, atype in zip(molFF.atoms, atypes):
        if(atom.atomic_number < 1):
            raise RuntimeError("non-real atoms are not supported")
        acharges.append(atom.charge)

        if atom.id not in keep_ids: continue
        a = res.addAtom()
        idmap[atom.id] = a
        if atom.id in selected_ids:
            a.name = atom.name
            a.atomic_number = atom.atomic_number
            a.mass = atom.mass
            a.charge = atom.charge
            tsys.setTypes(a, atype, atype)
        else:
            a.name = "$"+str(iexternal)
            a.atomic_number = -1
            iexternal += 1
        pdict={"type" : atype, "amu" : atom.mass, "memo" : ""}
        viparrConversionHelper.addOrUpdateParameterData(viparrff, "mass", pdict, False, False)
    for bond in molFF.bonds:
        id0 = bond.first.id
        id1 = bond.second.id
        if id0 not in keep_ids or id1 not in keep_ids: continue
        idmap[id0].addBond(idmap[id1])
    viparrConversionHelper.addTemplateData(viparrff.typer, tsys, symQ)
    # import IPython
    # IPython.embed()


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
                    ids = [ a.id for a in atoms ]
                    if len(set(ids) & keep_ids) == 4:
                        tsys.addImproper([ idmap[aid] for aid in ids])
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
            for t in table.terms:
                if(len(t.atoms)!=2):
                    raise RuntimeError("Pair Terms are malformed: %d != 2"%(len(t.atoms)))
                qSigEpsRef = []

                aids = [a.id for a in t.atoms]
                if aids[0] not in selected_ids or aids[1] not in selected_ids: continue
                for aid in aids:
                    termRef = nbtable.term(aid)
                    qSigEpsRef.append((acharges[aid], termRef["sigma"], termRef["epsilon"]))
                qij, aij, bij = t.param["qij"], t.param["aij"], t.param["bij"]

                qijRef = qSigEpsRef[0][0]*qSigEpsRef[1][0]
                qComb0 = my_es_scales[-1]*qijRef
                if(not numpy.isclose(qij, qComb0, atol=1E-4)):
                    delta = qij - qComb0
                    scale = qij / qijRef
                    print(f"Pair Term for ids {aids} qij != esScale*qi*qj {qij} {qComb0}: (delta, apparentScale {delta} {scale}")
                sigComb = pow(0.5 * (qSigEpsRef[0][1] + qSigEpsRef[1][1]), 6)
                epsComb = 4.0 * my_vdw_scales[-1] * (qSigEpsRef[0][2] * qSigEpsRef[1][2])**0.5
                bComb = epsComb * sigComb
                aComb = bComb * sigComb
                if(not numpy.isclose(aij, aComb)):
                    delta = aij - aComb
                    scale = aij / (aComb / my_vdw_scales[-1])
                    print(f"Pair Term for ids {aids} aij != vdwScale*eps_ij*sigma_ij^12 {aij} {aComb}: (delta, apparentScale {delta} {scale}")
                if(not numpy.isclose(bij, bComb)):
                    delta = bij - bComb
                    scale = bij / (bComb / my_vdw_scales[-1])
                    print(f"Pair Term for ids {aids} bij != vdwScale*eps_ij*sigma_ij^12 {bij} {bComb}: (delta, apparentScale {delta} {scale}")
        elif(table.name.startswith("constraint_")):
            pass
        else:
            raise RuntimeError("Unsupported table for conversion: %s"%(table.name))
    if( (len(non_pseudo_bonds) != termCounts["stretch_harm"]) or
        (len(angles) != termCounts["angle_harm"]) or
        (len(dihedrals) != termCounts["dihedral_trig"]) ):
        raise RuntimeError("Consistency check failed. There are missing valence terms")

    return viparrff
