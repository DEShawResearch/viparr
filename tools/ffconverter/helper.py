from __future__ import print_function
import numpy as np
import msys, viparr

def initializeSimpleForcefield(es_excl, vdw_excl, nbfix=""):
    if(len(vdw_excl)!=3 or len(es_excl)!=3):
       raise RuntimeError("initializeSimpleForcefield requires vdw and es exclusion lists to be length 3.")
    typer=viparr.TemplateTyper()
    rules = viparr.Rules()
    rules.info = []
    rules.vdw_func = 'LJ12_6_sig_epsilon'
    rules.vdw_comb_rule = "ARITHMETIC/GEOMETRIC"
    rules.plugins=["exclusions"]
    rules.setExclusions(4, es_excl, vdw_excl)
    rules.nbfix_identifier=nbfix
    return viparr.Forcefield(rules, typer)

def createParamTableFromSchema(viparrff, tname):
    schemas={
        "stretch_harm": (("r0", "fc"), "bonds"),
        "ureybradley_harm": (("r0", "fc"), "ureybradley"),
        "angle_harm": (("theta0", "fc"), "angles"),
        "dihedral_trig": (("phi0","fc0","fc1","fc2","fc3","fc4","fc5","fc6"),"propers"),
        "improper_trig": (("phi0","fc0","fc1","fc2","fc3","fc4","fc5","fc6"),"impropers"),
        "improper_harm": (("phi0","fc"),"impropers"),
        "improper_anharm": (("fc2","fc4"),"impropers"),
        "mass": (("amu",), "mass"),
        "vdw1": (("sigma","epsilon"),"vdw1"),
        "vdw1_14": (("sigma","epsilon"),"vdw1"),
        "vdw2": (("sigma","epsilon"), "vdw2"),
        "torsiontorsion_cmap": (("cmapid",),"cmap"),
        "virtuals_lc2": (("c1",),"virtuals"),
        "virtuals_lc2n": (("c1",),"virtuals")
    }
    pnames,plugin=schemas[tname]
    # make sure the plugin is in the current forcefield
    plugins= viparrff.rules.plugins
    if(plugin not in plugins):
        plugins.append(plugin)
        viparrff.rules.plugins=plugins

    # if the parameter table already exists, we are done
    if(viparr.Forcefield.HasParamTable(tname)): return

    # create and initialize the parameter table
    pt = msys.CreateParamTable()
    viparr.Forcefield.AddParamTable(tname, pt)
    pt.addProp('type', str)
    for n in pnames:
        if(n=="cmapid"):
            pt.addProp(n, int)
        else:
            pt.addProp(n, float)
    pt.addProp('memo', str)
    if("vdw" in tname):
        pt.addProp("nbfix_identifier", str)


def checkParameterValues(param, pdict):
    s=""
    keys=set(pdict.keys())-set(["memo","type"])
    for k in keys:
        v=pdict[k]
        if(param[k] != v):
           s+="   %s (%s -> %s)\n"%(k,str(param[k]),str(v))
    return s


def addOrUpdateParameterData(viparrff, tname, pdicts, modify, multiterm, checkReversed=True):
    createParamTableFromSchema(viparrff,tname)

    # multiterm parameters must come in as a group and have len>1
    assert( (isinstance(pdicts, list) and multiterm == True and len(pdicts)>0 and isinstance(pdicts[0], dict)) or
            (isinstance(pdicts, dict) and multiterm == False) )
    assert(not (modify and multiterm))

    if(isinstance(pdicts, list)):
        atypes=pdicts[0]["type"]
        for p in pdicts: assert(atypes==p["type"])
    else:
        atypes=pdicts["type"]

    if(checkReversed):
        atypes2=" ".join(reversed(atypes.split()))
        if(atypes2 < atypes): atypes, atypes2 = atypes2, atypes

        found = viparrff.findParams(tname, type=atypes)
        found2 = viparrff.findParams(tname, type=atypes2) if(atypes2 != atypes) else []
        if( len(found)>=0 and len(found2)==0 ):
            pass # atypes==atypes
        elif( len(found1) > 0 and len(found2) >0):
            err="Original and reversed parameter types both found in forcefield: "
            err+="table %s  new_entry %s  found_entry %s  found entry %s"%(tname, str(pdicts), str(found), str(found2))
            raise RuntimeError(err)
        else: # len(found2)>0 and len(found)==0
            err="Warning: Reversed parameter found in forcefield file... Continuing to use it: "
            err+="table %s  new_entry %s  found_entry %s"%(tname, str(pdicts), str(found2))
            print(err)
            atypes=atypes2
            found=found2
    else:
        found=viparrff.findParams(tname, type=atypes)

    if(modify):
        pdicts["type"]=atypes
        if(len(found)>1):
            raise UserWarning("Must have zero or one parameter in table %s with given key %s in order to modify: %s %s"%(tname,atypes,str(found),str(pdicts)))
        elif(len(found)==1):
            param=found[0]
        else:
            param = viparr.Forcefield.ParamTable(tname).addParam()
            viparrff.appendParam(tname,param)
        pdicts["type"]=atypes
        for k, v in pdicts.items(): param[k]=v
    else:
        if(isinstance(pdicts, dict)): pdicts=[pdicts]
        for p in pdicts: p["type"]=atypes

        if(len(found)==0):
            # add new set of parameters for the group
            for pdict in pdicts:
                param = viparr.Forcefield.ParamTable(tname).addParam()
                viparrff.appendParam(tname,param)
                for k, v in pdict.items(): param[k]=v
        else:
            assert(len(found)==len(pdicts))
            if(multiterm):
                ids=[p.id-found[0].id for p in found]
                if(ids != list(range(len(ids)))):
                    raise UserWarning("Non-contiguous multiterm parameter (duplicated key) detected for table %s with atypes %s: %s"%(tname,atypes,str(found)))
            s=""
            for f, pdict in zip(found, pdicts):
                s+=checkParameterValues(f, pdict)
            if(s!=""):
                s="Duplicate key with different values detected for table %s with atypes %s: \n%s"%(tname, atypes, s)
                raise UserWarning( s )
            return


def addVirtualSiteData(viparrff, vtype, pdict):
    tname="virtuals_"+vtype
    if(not viparr.Forcefield.HasParamTable(tname)): createParamTableFromSchema(viparrff,tname)

    atypes=pdict["type"]
    currentParams=viparrff.params(tname)
    count=0
    for p in viparrff.params(tname):
        ptypes=p["type"].split()
        ptypes,pset=" ".join(ptypes[:-1]),ptypes[-1]
        if(ptypes!=atypes): continue
        assert("pset%d"%(count)==pset)
        count+=1
        keys=set(p.keys())-set(["memo","type"])
        found=True
        for k in keys:
            found&= p[k] == pdict[k]
        if(found): return pset,p
    param = viparr.Forcefield.ParamTable(tname).addParam()
    viparrff.appendParam(tname,param)
    for k,v in pdict.items():
        param[k]=v
    pset="pset%d"%(count)
    param["type"]=atypes+" "+pset
    return pset,param

def addTemplateData(typer, newTemplate, symQ):

    res=newTemplate.system
    matches,err=typer.matchFragment(newTemplate,res.atoms)
    if(len(matches)>0):
        print("Duplicate template(s) found. Skipping current", newTemplate.name, [t[0].name for t in matches])
        return

    if(res.natoms>1 and symQ):
        atomic_charges=np.array([a.charge for a in res.atoms])
        ndigit=[]
        for q in atomic_charges:
            for i in xrange(10):
               if(round(q,i)==q):
                   ndigit.append(i)
                   break
            else:
               ndigit.append(i)
        ndigit=max(ndigit)

        formal_charges=np.full(res.natoms, sum(atomic_charges)/res.natoms)
        fc, bcis, bondChargeToAtomCharge=getFragmentPartitionedCharges(res, formal_charges, atomic_charges, None)
        cleanQ=fc+np.dot(bondChargeToAtomCharge,bcis)
        for atom, charge in zip(res.atoms, cleanQ):
            atom.charge = round(charge, ndigit+1)
    typer.addTemplate(newTemplate)


def getTopologyMatrix(mol, atypes=None):
    from collections import defaultdict
    assert(sorted([a.id for a in mol.atoms]) == list(range(mol.natoms)))
    assert(sorted([b.id for b in mol.bonds]) == list(range(mol.nbonds)))

    if(atypes is None):
        atypes = msys.ComputeTopologicalIds(mol)
    assert(len(atypes) == mol.natoms)
    # this works regardless of what type of data atypes is
    tidToIdx = {}
    tidToAtoms = defaultdict(list)
    for atom in mol.atoms:
        tid = atypes[atom.id]
        idx = len(tidToIdx)
        idx = tidToIdx.setdefault(tid, idx)
        tidToAtoms[tid].append(atom)

    nUniqueTids = len(tidToIdx)
    topology = np.zeros((mol.natoms, nUniqueTids))
    topologyInverse = np.zeros((nUniqueTids, mol.natoms))
    for tid, idx in tidToIdx.items():
        atoms = tidToAtoms[tid]
        factor = 1.0/float(len(atoms))
        for atom in tidToAtoms[tid]:
            topology[atom.id, idx] = 1.0
            topologyInverse[idx, atom.id] = factor

    bondChargeToAtomCharge = np.zeros((mol.natoms, mol.nbonds))
    for bond in mol.bonds:
        a0, a1 = bond.first, bond.second
        # positive bcis for atom with smaller atomtype
        sign=1.0 if(atypes[a0.id] < atypes[a1.id]) else -1.0
        bondChargeToAtomCharge[a0.id,bond.id] += sign
        bondChargeToAtomCharge[a1.id,bond.id] -= sign

    return topology, topologyInverse, bondChargeToAtomCharge


# split atomic partial charges q and atomic formal charges fc (with sum(q)=sum(fc)) into bcis
#      q=fc+T*bci
#    bci=T^-1(q-fc)
#  We can enforce topologically equivalent atomic charges via the transform:
#     P*q=P(fc+T*bci)
#     bci=(P*T)^-1*P*(q-fc)
#  For molecules, sum(formal_charges) == integer
def getMoleculePartitionedCharges(mol, formal_charges=None, atypes=None):
    nreal=len([1 for a in mol.atoms if a.atomic_number>0])

    atomChargesOrig = np.array([a.charge for a in mol.atoms])
    sumAC = sum(atomChargesOrig)
    sumFC = int(round(sumAC,0))
    excess = sumAC - sumFC
    if(np.abs(excess) > 1E-3):
        s="Molecule %s has non-integer total charge (%f). Partitioning will NOT be performed\n"%(mol.name,sumAC)
        s+="  Note: partitioning can be performed by calling getFragmentPartitionedCharges for fragments"
        raise RuntimeError(s)

    qRef = excess / nreal
    if(np.abs(qRef) > np.finfo(float).eps):
        if(np.abs(qRef) > 1E-5):
            print("Warning: integerization of atomic charges for '%s' results in large perturbation per atom."%(mol.name))
            print("         '%s' excess per atom: %e ( %e / %d )"%(mol.name, qRef, excess, nreal))
        #modify atom charges so that fragment has integer charge
        atomChargesOrig -= np.array([qRef if a.atomic_number>0 else 0.0 for a in mol.atoms  ])
    sumAC = float(sumFC)

    if(formal_charges is None):
        formalChargesOrig = np.array([a.formal_charge for a in mol.atoms], dtype=float)
    else:
        formalChargesOrig = np.array(formal_charges,dtype=float)
    if(sum(formalChargesOrig) != sumFC):
        fakeAtoms = [a for a in mol.atoms if a.atomic_number<1]
        for a in fakeAtoms: a.atomic_number=1
        msys.AssignBondOrderAndFormalCharge(mol,sumFC)
        formalChargesOrig=np.array([a.formal_charge for a in mol.atoms], dtype=float)
        for a in fakeAtoms: a.atomic_number=-1

    return getFragmentPartitionedCharges(mol, formalChargesOrig, atomChargesOrig, atypes)


def getFragmentPartitionedCharges(mol, formal_charges, atomic_charges, atypes):
    """ robustly split atomic partial charges "q" into atomic formal charges "fc"
    (with sum(q) == sum(fc)) and bond charge increments "bci"s, thus:
        q = fc + T*bci
        bci = T^-1(q-fc)
    where T [bondChargeToAtomCharge] is a mapping matrix that takes bcis to atoms
    We can enforce topologically equivalent atomic charges (and bcis) via the transform:
        P*q = P(fc + T*bci)
        bci = (P*T)^-1*P*(q-fc)
    where P [topologyInverse] maps atomic charges to topologically unique charges

    This function does not require that sum(formal_charges) is integer, the only
    consistency check is that sum(formal_charges) == sum(atomic_charges)
    Args:
        mol (msys system)
        formal_charges: numpy array (len natom) of the atomic formal charges
        atomic_charges: numpy array (len natom) of the atomic charges.
        atypes: list (len natom) of the type for each atom to determine equivalence

    Returns:
        fc: numpy array (len natom) of the rounded and topologically idential formal charges
        bci: numpy array (len bonds) of the rounded bond charge increments
        bondChargeToAtomCharge: matrix (natom x nbond) to transform from bcis to atomic charges
        """

    tol=1E-5
    assert(abs(sum(formal_charges) - sum(atomic_charges)) <= tol)

    topology, topologyInverse, bondChargeToAtomCharge = getTopologyMatrix(mol, atypes)

    #convert to topologically unique and rounded atomic charge vectors:
    formalUniqueCharges = np.dot(topologyInverse, formal_charges)
    atomicUniqueCharges = np.dot(topologyInverse, atomic_charges)

    bondToUniqueAtom = np.dot(topologyInverse, bondChargeToAtomCharge)

    bcis, _, rank, s = np.linalg.lstsq(bondToUniqueAtom, atomicUniqueCharges-formalUniqueCharges, rcond=1E-8)
    #print "BCI: ",mol,bcis
    #print "RANK,S:  ",mol,rank,len(atomUniqueCharges)-1,s

    fc = np.dot(topology, formalUniqueCharges)
    q = fc + np.dot(bondChargeToAtomCharge, bcis)
    # make sure our bci based solution reproduces the symmetry equivalent one
    atomCharges = np.dot(topology, atomicUniqueCharges)
    if( (np.abs(q-atomCharges) > tol).any() ):
        print("Decomposed charges and reference charges are not within tolerance")
        for a,b in zip(q,atomCharges):
            print(a, b, a-b)
        raise RuntimeError("Bad charge decomposition")

    return fc, bcis, bondChargeToAtomCharge

