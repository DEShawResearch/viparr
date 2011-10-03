from __future__ import print_function

import os
import sys
from argparse import Namespace
from collections import OrderedDict,defaultdict
import numpy as np
import msys,viparr
from viparr import ffconverter

convaname={
  "DISU" : {"1CB":"CB","1SG":"SG","1HG1":"HG1","2CB":"BAD1","2SG":"BAD2","2HG1":"BAD3" },
  "TP1"  : {"1HH":"HH"},
  "TP2"  : {"1HH":"HH"},
  "SP1"  : {"1HG1":"HG1"},
  "SP2"  : {"1HG1":"HG1"},
  "THP1" : {"1HG1":"HG1"},
  "THP2" : {"1HG1":"HG1"},
  "Zn2+" : {"ZN" : "Zn2+"},
  "Cs+"  : {"CES" : "Cs+"},
  "Mg2+" : {"MG"  : "Mg2+"},
  "Cl-"  : {"CLA" : "Cl-"},
  "Ca2+" : {"CAL" : "Ca2+"},
  "Cd2+" : {"CD" : "Cd2+"},
  "Ba2+" : {"BAR" : "Ba2+"},
  "Li+"  : {"LIT" : "Li+"},
  "Na+"  : {"SOD" : "Na+"},
  "K+"   : {"POT" : "K+"},
  "Rb+"  : {"RUB" : "Rb+"}
  }

# Better residue names
convres={
  "HSE" : "HIE",
  "HSP" : "HIP",
  "HSD" : "HIS",
  "ZN2" : "Zn2+",
  "MG"  : "Mg2+",
  "CD2" : "Cd2+",
  "CLA" : "Cl-",
  "CAL" : "Ca2+",
  "BAR" : "Ba2+",
  "LIT" : "Li+",
  "SOD" : "Na+",
  "POT" : "K+",
  "RUB" : "Rb+",
  "CES" : "Cs+"
  }

def clean_line(l):
    l=l.strip()
    cstart=l.find('!')
    comment=""
    if(cstart<0):
        words=l.split()
    else:
        comment=l[cstart+1:].strip()
        words=l[:cstart].split()

    # viparr1 uses python to eval the json data.
    # Python comments in eval'd strings cause problems
    comment=comment.replace('#','_')
    return(words,comment)

def sec_parse_from_re(lines,regex):
    parsed=[]
    for l in lines:
        lclean=l.strip()
        if(len(lclean)==0): break
        match=regex.match(lclean)
        if(not match):
            print("Failed to match regex to line: ", lclean)
            print("  Regex=",regex.pattern)
            assert(False)
        parsed.append(["" if m is None else m for m in match.groups()])
    return parsed

def findResidueAtom(res,atomName):
    atoms=[ a for a in res.atoms if a.name==atomName]
    assert(len(atoms)<=1)
    return None if len(atoms)==0 else atoms[0]

def findImproper(tsys,improperName):
    res=tsys.system.residue(0)
    improper=[findResidueAtom(res,aname) for aname in improperName]
    if None in improper: return None
    impropers=[imp for imp in tsys.impropers if imp==improper]
    assert(len(impropers)<=1)
    return None if len(impropers)==0 else impropers[0]

def removeOpAtoms(operations,removeNames):
    resname,res_is_patch,resq=operations["resprops"]

    operations["addAtoms"]=[ entry for entry in operations["addAtoms"] if entry[0] not in removeNames ]
    operations["deleteAtoms"]=[ aname for aname in operations["deleteAtoms"] if aname not in removeNames ]
    kill=set(removeNames)
    operations["addBonds"]=[ entry for entry in operations["addBonds"] if len(set(entry)&kill)==0 ]
    operations["deleteBonds"]=[ entry for entry in operations["deleteImpropers"] if len(set(entry)&kill)==0 ]
    operations["addImpropers"]=[ entry for entry in operations["addImpropers"] if len(set(entry)&kill)==0 ]
    operations["deleteImpropers"]=[ entry for entry in operations["deleteImpropers"] if len(set(entry)&kill)==0 ]
    operations["addCmaps"]=[ entry for entry in operations["addCmaps"] if len(set(entry)&kill)==0 ]
    operations["deleteCmaps"]=[ entry for entry in operations["deleteCmaps"] if len(set(entry)&kill)==0 ]

    resq=sum([entry[2] for entry in operations["addAtoms"]])
    operations["resprops"]=[resname,res_is_patch,resq]

def findOpAtomIndex(operations,atomName):
    idx=[i for i,a in enumerate(operations["addAtoms"]) if a[0]==atomName]
    assert(len(idx)<=1)
    return None if len(idx)==0 else idx[0]

def fix_aname(an,res):
    global convaname

    if (isinstance(an,str)):
        names=[an,]
    else:
        names=an

    newnames={}
    if res in convaname:
        newnames=convaname[res]
    for i in range(len(names)):
        names[i]=names[i].upper()
        if names[i] in newnames:
            names[i]=newnames[names[i]]

    if (isinstance(an,str)):
        an=names[0]
    else:
        an=names
    return an


def preload_merge(viparrdir,FFconv):
    if not os.path.exists(viparrdir) or not os.path.isdir(viparrdir):
        raise UserWarning("Bad base directory for viparr merge: "+viparrdir)
    file=os.path.join(viparrdir,"mass")
    if not os.path.exists(file):
        raise UserWarning("Couldnt find mass table: "+file)
    fh=open(file,"r")
    lines=fh.read()
    fh.close()
    data=eval(lines)
    if(len(data)==0): return
    for d in data:
        if isinstance(d,list):
            atype=d[0]
            mass=d[1]
        elif isinstance(d,dict):
            atype=d["type"]
            mass=d["params"]["amu"]
        else:
            raise UserWarning("Bad mass file: %s \n  contents: %s"%(file,d))

        FFconv.atype2atomicnumber[atype]=msys.GuessAtomicNumber(mass)


def parse_stretch(FFconv,lines):
    for l in lines[1:]:
        (words,comment)=clean_line(l)
        if(len(words)==0): continue
        params={
            "type": " ".join(FFconv.fix_atypes(words[0:2])),
            "fc"  : float(words[2]),
            "r0"  : float(words[3]),
            "memo": comment.strip()
        }
        ffconverter.addOrUpdateParameterData(FFconv.viparrff,"stretch_harm", params, False, False)


def parse_angle(FFconv,lines):
    for l in lines[1:]:
        (words,comment)=clean_line(l)
        if(len(words)==0): continue
        atypes=" ".join(FFconv.fix_atypes(words[0:3]))
        params={
            "type"  : atypes,
            "fc"    : float(words[3]),
            "theta0": float(words[4]),
            "memo"  : comment.strip()
        }
        ffconverter.addOrUpdateParameterData(FFconv.viparrff,"angle_harm", params, False, False)

        if(len(words) == 7 ):
            params={
                "type"  : atypes,
                "fc"    : float(words[5]),
                "r0"    : float(words[6]),
                "memo"  : comment.strip()
            }
            ffconverter.addOrUpdateParameterData(FFconv.viparrff,"ureybradley_harm", params, False, False)

def parse_proper(FFconv,lines):

    # combine all the torsion fc's by type, then add grouped parameters to forcefield
    torsionList=[]
    for l in lines[1:]:
        (words,comment)=clean_line(l)
        if(len(words)==0): continue

        for i in range(4):
            words[i]=words[i].replace('X','*')

        fc     = float(words[4])
        pn     = int(words[5])
        multi  ="fc%d"%(pn)
        phase  = float(words[6])
        atypes=" ".join(FFconv.fix_atypes(words[0:4]))

        possible=[]
        for params in reversed(torsionList):
            if(params["type"]!=atypes): break
            possible.append(params)
        possible=possible[::-1]

        found=False
        for params in possible:
            if( (phase == 0.0 or phase == 180.0) and params["phi0"] == 0.0 ):
                found=True
                if params[multi] != 0.0:
                    #if( (phase == 0.0 and data[1+pn]==-fc) or (phase == 180.0 and data[1+pn]==-fc)): return
                    raise UserWarning("multiplicity has already been assigned a value: %s %d %f %f"%(str(atypes),pn,params[multi],fc))
                params["fc0"]=float("%.6f"%(params["fc0"]+fc))
                if (phase==0.0):
                    params[multi]=fc
                elif (phase==180.0):
                    params[multi]=0.0 if 0.0==fc else -fc
                if comment not in params["memo"]: params["memo"]=params["memo"]+" "+comment
                break
            # merge for phase != {0,180}
            elif(phase == params["phi0"]):
                found=True
                if params[multi] != 0.0:
                    #if(data[1+pn]==fc): return
                    raise UserWarning("multiplicity has already been assigned a value: %s %d %f %f"%(str(atypes),pn,params[multi],fc))
                params["fc0"]+=fc
                params[multi]=fc
                if comment not in params["memo"]: params["memo"]=params["memo"]+" "+comment
                break
        if(found): continue

        params={"type":atypes, "phi0":0.0,
                "fc0":0.0, "fc1":0.0, "fc2":0.0,
                "fc3":0.0, "fc4":0.0, "fc5":0.0,
                "fc6":0.0, "memo":comment}
        params["fc0"]=fc
        if (phase==0.0):
            params[multi]=fc
        elif (phase==180.0):
            params[multi]=0.0 if 0.0==fc else -fc
        else:
            params["phi0"]=phase
            params[multi]=fc
        torsionList.append(params)

    groupedList=[]
    for params in torsionList:
        if(len(groupedList) and groupedList[-1][0]["type"]==params["type"]):
            groupedList[-1].append(params)
        else:
            groupedList.append([params])

    for params in groupedList:
        ffconverter.addOrUpdateParameterData(FFconv.viparrff,"dihedral_trig", params, False, True)

def parse_improper(FFconv,lines):
    for l in lines[1:]:
        (words,comment)=clean_line(l)
        if(len(words)==0): continue

        for i in range(4):
            words[i]=words[i].replace('X','*')
        params={
            "type"  : " ".join(FFconv.fix_atypes(words[0:4])),
            "fc"    : float(words[4]),
            "phi0"  : float(words[6]),
            "memo"  : comment.strip()
        }
        ffconverter.addOrUpdateParameterData(FFconv.viparrff,"improper_harm", params, False, False, False)

def parse_cmap(FFconv,lines):
    if(FFconv.remove_cmap): return

    cmapNtotal=0
    cmapTable=None
    cmapPhi=None
    cmapPsi=None
    for l in lines[1:]:
        (words,comment)=clean_line(l)
        if(len(words)==0): continue

        if(cmapNtotal==0):
            if(len(words) !=9 ):
                raise UserWarning("Can't parse cmap term: "+str(words))

            t1,t2=words[0:4],words[4:8]
            id=len(FFconv.viparrff.cmap_tables)+1
            params={
                "type":" ".join(t1+t2),
                "cmapid":id,
                "memo":comment
            }
            ffconverter.addOrUpdateParameterData(FFconv.viparrff,"torsiontorsion_cmap",params,False, False, False)

            cmapTable = msys.CreateParamTable()
            cmapTable.addProp('phi', float)
            cmapTable.addProp('psi', float)
            cmapTable.addProp('energy', float)
            FFconv.viparrff.addCmapTable(cmapTable)

            npts=int(words[8])
            assert(npts==24)
            spacing=15 #360.0/npts
            cmapNtotal = npts*npts
            cmapPhi = -180.0
            cmapPsi = -180.0
            continue

        for t in map(float,words):
            param = cmapTable.addParam()
            for k,v in [("phi", cmapPhi),("psi", cmapPsi),("energy", t)]:
                param[k]=v

            cmapPsi += spacing
            if(cmapPsi > 180.0-0.5*spacing):
                cmapPsi = -180.0
                cmapPhi += spacing
            cmapNtotal-=1
            assert(cmapNtotal>=0)

def parse_lj(FFconv,lines):
    for l in lines[1:]:

        (words,comment)=clean_line(l)
        if(len(words)==0): continue

        params={
            "type" : " ".join(FFconv.fix_atypes(words[0:1])),
            "sigma":  float(words[3])*FFconv.conv_lj[0],
            "epsilon": float(words[2])*FFconv.conv_lj[1],
            "nbfix_identifier": FFconv.nbfixID,
            "memo" : comment
        }
        ffconverter.addOrUpdateParameterData(FFconv.viparrff,"vdw1", params, False, False)

        if(len(words) == 7 ):
            #del(params["nbfix_identifier"])
            params["sigma"]= float(words[6])*FFconv.conv_lj[0]
            params["epsilon"]=float(words[5])*FFconv.conv_lj[1]
            ffconverter.addOrUpdateParameterData(FFconv.viparrff,"vdw1_14", params, False, False)

def parse_nop(FFconv,lines):
    return

def parse_NBFIX(FFconv,lines):
    for l in lines[1:]:
        words,comment=clean_line(l)
        nwords=len(words)
        if nwords==0: continue
        if(nwords!=4):
            print("Invalid NBFIX statement (nwords=%d): %s" %(nwords,str(words)))

        params={
            "type" : " ".join(FFconv.fix_atypes(words[0:2])),
            "sigma":  float(words[3])*FFconv.conv_nbfix_lj[0],
            "epsilon": float(words[2])*FFconv.conv_nbfix_lj[1],
            "nbfix_identifier": FFconv.nbfixID,
            "memo" : comment
        }
        ffconverter.addOrUpdateParameterData(FFconv.viparrff,"vdw2", params, False, False)


def parse_mass(FFconv,lines):
    for l in lines:
        (words,comment)=clean_line(lines[0])
        nwords=len(words)
        if(nwords==0): continue
        if(nwords!=4 and nwords !=5):
            print("Invalid MASS statement (nwords=%d): %s"%(nwords,str(words)))
            return

        idx=int(words[1])
        type=FFconv.fix_atypes(words[2])
        mass=float(words[3])
        if(mass==0.0):
            FFconv.atype2atomicnumber[type]=0
        else:
            FFconv.atype2atomicnumber[type]=msys.GuessAtomicNumber(mass)

        params={
            "type": type,
            "amu": mass,
            "memo": comment
        }
        ffconverter.addOrUpdateParameterData(FFconv.viparrff,"mass", params, False, False)


def build_residue(FFconv,operations,_tsys=None):
    amap={}
    tsysCharge=0.0

    initialName=None
    if(_tsys is None):
        assert("deleteAtoms" not in operations and "deleteImpropers" not in operations)
        tsys = viparr.TemplatedSystem()
        res = tsys.system.addResidue()
    else:
        tsys=_tsys.clone()
        res=tsys.system.residue(0)
        initialName=res.name
        for atomName in operations["deleteAtoms"]:
            atom=findResidueAtom(res,atomName)
            if(atom is None):
                print("couldnt find atom to remove: ",initialName,atomName)
                continue
            atom.remove()
            for improper in tsys.impropers:
                if(atom in improper): tsys.removeImproper(improper)
            for cmap in tsys.cmaps:
                if(atom in cmap): tsys.removeCmaps(cmap)
        for improperName in operations["deleteImpropers"]:
            improper=findImproper(tsys,improperName)
            if(improper is None): continue
            tsys.removeImproper(improper)
        for a in res.atoms:
            tsysCharge+=a.charge
            amap[a.name]=a

    resname,res_is_patch,resq=operations["resprops"]
    aqsum=0.0
    res.name=resname
    for aname, anum, acharge, atype, pset in operations["addAtoms"]:
        a=amap.get(aname,None)
        if(a is None):
            a=res.addAtom()
            amap[aname]=a
        else:
            tsysCharge-=a.charge
        a.name=aname
        a.atomic_number=anum
        a.charge=acharge
        aqsum+=acharge
        if(anum >= 0): tsys.setTypes(a,atype,atype,pset)

    for anames in operations["addBonds"]:
        aname0,aname1=[FFconv.bsplit.get(n,n) for n in anames]
        assert(not (aname0[0] == '$' and aname1[0] == '$'))

        if aname0[0] =='$' and aname1[0] == '$':
            raise UserWarning("Cant have bond that isnt part of this residue")
        elif aname0[0] =='$':
            assert(aname0 not in amap)
            a=res.addAtom()
            a.name=aname0
            a.atomic_number=-1
            amap[a.name]=a
        elif aname1[0] =='$':
            assert(aname1 not in amap)
            a=res.addAtom()
            a.name=aname1
            a.atomic_number=-1
            amap[a.name]=a

        a0=amap[aname0]
        a1=amap[aname1]
        a0.addBond(a1)

    for improperNames in operations["addImpropers"]:
        try:
            imprAtoms=[amap[FFconv.bsplit.get(n,n)] for n in improperNames]
            tsys.addImproper(imprAtoms)
        except:
            print("Couldnt add improper %s to residue %s"%(str(improperNames),resname))

    for cmapNames in operations["addCmaps"]:
        cmapAtoms=[amap[FFconv.bsplit.get(n,n)] for n in cmapNames]
        tsys.addCmap(cmapAtoms)

    for vtype, anames in operations["defineVirtuals"]:
        vname="virtual_"+vtype
        tsys.addPseudoType(vname, len(anames))
        virtAtoms=[amap[FFconv.bsplit.get(n,n)] for n in anames]
        virtAtoms[0].addBond(virtAtoms[1])
        tsys.addPseudoSites(vname, virtAtoms)

    qBuildSum=tsysCharge+resq
    systemQ=sum([a.charge for a in res.atoms])
    assert np.isclose(qBuildSum,systemQ)
    if(res.name not in ['fe1a']):
        intQ=round(qBuildSum,0)
        if(np.abs(qBuildSum-intQ)>1E-5):
            print("NON-INTEGER template charge: ",tsys.name,initialName,res.name,qBuildSum,systemQ,intQ)
            #assert(False)
    fragments=tsys.system.updateFragids()
    if(len(fragments)>1):
        print("WARNING: Template atoms are disconnected",operations["resprops"])
        for i,frag in enumerate(fragments):
            print(i,[a.name for a in frag])
        #assert(False)
    return tsys


def parse_residue(FFconv, lines):
    from collections import Counter
    from itertools import chain

    otherres = ["+", "-"]
    operations = defaultdict(list)
    atomTypeMap = {}
    # unfortunatly, charmm only defines missing bonds in one direction.
    # we are left to try and infer when we need to add an external link 
    # Sometimes we add too many links. This attempts to generically find 
    # and fix those cases
    extraBonds = []
    for l in lines:
        (words, comment) = clean_line(l)
        if(len(words) == 0): continue

        if(words[0][0:4].upper() in ["IC", "ANGL","THET","DONO","ACCE", "BILD","PATC", "DIHE"]):
            continue
        elif(words[0][0:4].upper() == "RESI" or words[0][0:4].upper() == "PRES"):
            resname = words[1]
            if(resname in FFconv.skipres or (len(FFconv.onlyThese) and resname not in FFconv.onlyThese)): return
            res_is_patch = True if words[0][0:4].upper() == "PRES" else False
            resname = convres.get(resname,resname)
            if(len(words) > 2):
                resq = float(words[2])
            else:
                print("Warning: No charge found for residue %s. assuming 0.0"%(resname))
                resq = 0.0
            if(resname == "DDAOP"): resq = 1.0
            operations["resprops"].extend([resname, res_is_patch,resq])
        elif(words[0][0:4].upper() == "GROU"):
            continue
        elif(words[0][0:4].upper() == "ATOM"):
            aname = fix_aname(words[1],resname)
            atype = FFconv.fix_atypes(words[2])
            acharge = float(words[3])
            anum = FFconv.atype2atomicnumber.get(atype, None)

            if(anum is None):
                print("Error while generating template for residue: "+resname)
                print("  Atom type '%s'->'%s' wasnt defined in a MASS statement"%(words[2],atype))
                print("  You are most likely trying to convert an incomplete forcefield.")
                print("  Try using '--merge' or additional '-p' options to this conversion utility")
                print("  Currently known atomtypes and atomic numbers:")
                for k,v in FFconv.atype2atomicnumber.items(): print("    %s: %d"%(k,v))
                sys.exit(1)
            operations["addAtoms"].append((aname, anum, acharge, atype, ""))
            atomTypeMap[aname] = atype

        elif(words[0][0:4].upper() == "BOND" or words[0][0:4].upper() == "DOUB" or words[0][0:4].upper() == "TRIP"):
            anames = fix_aname(words[1:], resname)
            assert(len(anames)%2 == 0)
            for aname0, aname1 in zip(anames[::2], anames[1::2]):
                if(res_is_patch):
                    assert(not ('-' in aname0 or '+' in aname0 or '-' in aname1 or '+' in aname1))
                if(aname1[0] in otherres): aname0, aname1 = aname1, aname0
                if(aname0[0] in otherres):
                    operations["addBonds"].append((aname0, aname1))
                    prefix, aname0 = aname0[0], aname0[1:]
                    aname1 = otherres[(otherres.index(prefix)+1)%2] + aname1
                    operations["addBonds"].append((aname0, aname1))
                    extraBonds.append((aname0, aname1))
                else:
                    operations["addBonds"].append((aname0, aname1))

        elif(words[0][0:4].upper() == "IMPR" or words[0][0:4].upper() == "IMPH"):
            anames = fix_aname(words[1:], resname)
            assert(len(anames)%4 == 0)
            for improperNames in zip(anames[::4], anames[1::4], anames[2::4], anames[3::4]):
                operations["addImpropers"].append(improperNames)

        elif(words[0][0:4].upper() == "CMAP"):
            anames = fix_aname(words[1:], resname)
            assert(len(anames) == 8)
            operations["addCmaps"].append(anames)

        elif(words[0][0:4].upper() == "LONE"):
            lptype = words[1][0:4].upper()
            vtype = None
            if(lptype == "COLI"):
                if(len(words) == 7):
                    widx = [5]
                elif(len(words) == 9):
                    widx = [5, 7]
                else:
                    raise UserWarning("I dont know how to parse this lonepair: "+str(words[1:]))
                anames = fix_aname(words[2:5], resname)
                params={"c1": 0.0,
                        "memo": " ".join(words)
                }
                final = []
                for w in widx:
                    p = -float(words[w+1])
                    if (words[w][0:4].upper() == "DIST" and p != 0.0):
                        final.append(("lc2n", p))
                    elif(words[w][0:4].upper() == "SCAL" and p != 0.0):
                        final.append(("lc2", p))
                if(len(final) != 1):
                    UserWarning("I dont know how to parse this lonepair: "+str(words[1:]))
                vtype,params["c1"] = final[0]
            elif(lptype == "RELA"):
                assert(len(words) == 12)
                anames = fix_aname(words[2:6], resname)
                vtype = "fdat3"
                params = {"type":" ".join(anames),
                          "memo": " ".join(words)
                }
                raise UserWarning("I dont know how to parse this lonepair: "+str(words[1:]))
            else:
                raise UserWarning("I dont know how to parse this lonepair: "+str(words[1:]))
            params["type"] = " ".join(map(atomTypeMap.get,anames[1:]))
            subtype,refparam = ffconverter.addVirtualSiteData(FFconv.viparrff, vtype, params)
            idx = findOpAtomIndex(operations,anames[0])
            assert(idx is not None)
            (aname, anum, acharge, atype, pset) = operations["addAtoms"][idx]
            if(pset != ""):
                keys = set(refparam.keys()) - set(["type","memo"])
                print(subtype, pset, [(refparam[k],params[k]) for k in keys])
                assert(subtype == pset and np.all([refparam[k] == params[k] for k in keys]))
            else:
                operations["addAtoms"][idx] = (aname, anum, acharge, atype, subtype)
                operations["defineVirtuals"].append((vtype, anames))

        elif(words[0][0:4].upper() == "DELE"):
            if not res_is_patch: raise UserWarning("DELETE directive inside of non-patch")
            if (words[1][0:4].upper() == "ATOM" and len(words) == 3):
                anames=fix_aname(words[2:], resname)
                operations["deleteAtoms"].extend(anames)
            elif(words[1][0:4].upper() == "IMPR" and len(words) == 6):
                anames=fix_aname(words[2:6], resname)
                operations["deleteImpropers"].append(anames)
            elif(words[1][0:4].upper() == "ANGL"):
                pass
            elif(words[1][0:4].upper() != "ACCE"):
                s="Unknown DELETE directive " + l.strip()
                raise UserWarning(s)
        else:
            raise UserWarning("I dont know what to do with: " + str(words))

    addBonds = operations["addBonds"]
    count = Counter(chain.from_iterable(addBonds))
    filterN = [ b for b in extraBonds if (b[0]=='N' and count['N'] == 4) ]
    assert len(filterN) < 2
    for b in filterN: addBonds.remove(b)
    operations["addBonds"] = list(set([tuple(sorted(k)) for k in addBonds]))

    qres = operations["resprops"][2]
    qatoms = sum([entry[2] for entry in operations["addAtoms"]])
    if(not np.isclose(qres, qatoms)):
        s="Residue (%s) charge (%f) and sum over atoms charge (%f) are not close: %f"%(resname, qres, qatoms, qres-qatoms) 
        raise RuntimeError(s)

    if(res_is_patch):
        assert(resname not in FFconv.patches)
        FFconv.patches[resname] = operations
    else:
        tsys = build_residue(FFconv,operations)
        ffconverter.addTemplateData(FFconv.viparrff.typer, tsys, True)

def apply_pres(FFconv, prefix, pname, resname, presname, addToTemplates=True):

    tsys=FFconv.viparrff.typer.findTemplate(resname)
    operations=FFconv.patches.get(pname,None)
    assert(len(tsys)<=1)
    if(len(tsys)==0 or operations is None): return None

    tsys=tsys[0].clone()
    res=tsys.system.residue(0)

    if(pname=="DISU"):
        if("+SG" not in FFconv.bsplit):
            FFconv.bsplit["+SG"]="$%d"%(len(FFconv.bsplit)+1)
        removeOpAtoms(operations,["BAD1", "BAD2", "BAD3"])
        operations["addBonds"].append(("SG", "+SG"))

    # clean up residue
    if(prefix=='N'):
        atom=findResidueAtom(res,FFconv.bsplit["-C"])
        if(atom is not None): atom.remove()
    elif(prefix=='C'):
        atom=findResidueAtom(res,FFconv.bsplit["+N"])
        if(atom is not None): atom.remove()
    elif(prefix=='5'):
        atom=findResidueAtom(res,FFconv.bsplit["-O3'"])
        if(atom is not None): atom.remove()
    elif(prefix=='3'):
        atom=findResidueAtom(res,FFconv.bsplit["+P"])
        if(atom is not None): atom.remove()

    patched=build_residue(FFconv,operations,tsys)

    res=patched.system.residue(0)
    res.name=presname

    # add new residue
    if(addToTemplates):
        matches,err=FFconv.viparrff.typer.matchFragment(patched,res.atoms)
        if(len(matches)>0):
            print("Duplicate template(s) found. Skipping current",patched.name,[t[0] for t in matches])
        else:
            FFconv.viparrff.typer.addTemplate(patched)
    return patched

def mergeResiduesViaPatch(FFconv, resname1, resname2, pname, rsuffix):
    tsys1 = FFconv.viparrff.typer.findTemplate(resname1)
    tsys2 = FFconv.viparrff.typer.findTemplate(resname2)

    operations = FFconv.patches.get(pname,None)

    if(len(tsys1) == 0 or len(tsys2) == 0 or operations is None): return None
    tsys1 = tsys1[0].clone()
    tsys2 = tsys2[0].clone()
    residue1 = tsys1.system.residue(0)
    residue2 = tsys2.system.residue(0)
    print("MergeResidueOperations: ",operations)

    if(len(rsuffix) != 2 ): raise UserWarning("bad rename_suffix")
    if(rsuffix != ("","")):
        if(rsuffix[0] != ""):
            for atom in residue1.atoms: atom.name += rsuffix[0]
        if(rsuffix[1] != ""):
            for atom in residue2.atoms: atom.name += rsuffix[1]
        entries = []
        for (aname, anum, acharge, atype, memo) in operations["addAtoms"]:
            if(aname[0] == "1"):
                aname = aname[1:] + rsuffix[0]
            elif(aname[0] == "2"):
                 aname = aname[1:] + rsuffix[1]
            entries.append((aname, anum, acharge, atype, memo))
        operations["addAtoms"] = entries
        entries = []
        for b in operations.get("addBonds",[]):
            bond = []
            for aname in b:
                if(aname[0] == "1"):
                    aname = aname[1:] + rsuffix[0]
                elif(aname[0] == "2"):
                    aname = aname[1:] + rsuffix[1]
                bond.append(aname)
            entries.append(tuple(bond))
        if(len(entries)): operations["addBonds"] = entries
        entries = []
        for aname in operations.get("deleteAtoms", []):
            if(aname[0] == "1"):
                aname = aname[1:] + rsuffix[0]
            elif(aname[0] == "2"):
                aname = aname[1:] + rsuffix[1]
            entries.append(aname)
        if(len(entries)): operations["deleteAtoms"] = entries
    amap = {}
    for atom in residue2.atoms:
        newAtom = residue1.addAtom()
        newAtom.atomic_number = atom.atomic_number
        newAtom.charge = atom.charge
        newAtom.name = atom.name
        tsys1.setTypes(newAtom, tsys2.btype(atom), tsys2.nbtype(atom), tsys2.pset(atom))
        amap[atom] = newAtom
    for bond in tsys2.system.bonds:
        amap[bond.first].addBond(amap[bond.second])

    patched = build_residue(FFconv, operations, tsys1).clone()

    return patched

def fixup_amino_acid_patches(FFconv):
    pconv={"ACE":"ACE","CT3":"NMA"}

    for pname,rname in pconv.items():
        if(pname not in FFconv.patches): continue
        operations =FFconv.patches[pname]

        if("deleteAtoms" in operations):
            raise UserWarning("Cant convert patch %s to residue %s with atoms to delete"%(pname,rname))

        # clear cmap (its handled by the main aa residues)
        if("addCmaps" in operations):
            del(operations["addCmaps"])

        if(rname=='ACE'):
            # This is handled by main aa residues
            operations["addBonds"].remove(("CY","N"))
            operations["addBonds"].append(("CY","+N"))
            operations["addImpropers"].remove(("N", "CY", "CA", "HN"))
            operations["addImpropers"].remove(('CY', 'CAY', 'N', 'OY'))
            operations["addImpropers"].append(('CY', 'CAY', '+N', 'OY'))

        if(rname=='NMA'):
            # This is handled by main aa residues
            operations["addAtoms"]=[ entry for entry in operations["addAtoms"] if entry[0] not in ["C","O"] ]
            operations["addBonds"]=[ entry for entry in operations["addBonds"] if "C" not in entry and "O" not in entry ]
            operations["addImpropers"]=[ entry for entry in operations["addImpropers"] if "C" not in entry and "O" not in entry ]
            operations["addBonds"].append(("NT",'-C'))
            operations["addImpropers"].append(("NT",'-C',"CAT","HNT"))

        tsys=build_residue(FFconv,operations)
        tsys.resname=rname
        del(FFconv.patches[pname])
        ffconverter.addTemplateData(FFconv.viparrff.typer, tsys, True)


def patch_amino_acids(FFconv):
    #This handles ACE/ACP (same patch) and CT3
    fixup_amino_acid_patches(FFconv)

    apply_pres(FFconv,"", "DISU", "CYS", "CYX" )

    # Protonation States (maestro format):
    # GLH = GLU + H  [GLUP]
    # ASH = ASP + H  [ASPP]
    # LYN = LYS - H  [LSN]
    # ARN = ARG - H  [????]

    apply_pres(FFconv,"", "ASPP", "ASP", "ASH")
    apply_pres(FFconv,"", "GLUP", "GLU", "GLH")
    apply_pres(FFconv,"", "LSN",  "LYS", "LYN" )
    apply_pres(FFconv,"", "SERD", "SER", "dSER" )
    # phosphorylated amino acid compounds
    apply_pres(FFconv,"", "TP1", "TYR", "MPTYR")
    apply_pres(FFconv,"", "TP2", "TYR", "DPTYR")
    apply_pres(FFconv,"", "SP1", "SER", "MPSER")
    apply_pres(FFconv,"", "SP2", "SER", "DPSER")
    apply_pres(FFconv,"", "THP1", "THR", "MPTHR")
    apply_pres(FFconv,"", "THP2", "THR", "DPTHR")


    #apply NTER and CTER pres to aa residues:
    aa=['ALA','ARG','ASN','ASP','CYS',
        'GLN','GLU','GLY','HIS','HIE',
        'HIP','ILE','LEU','LYS','MET',
        'PHE','PRO','SER','THR','TRP',
        'TYR','VAL','CYM',
        'ASH','GLH','LYN','dSER',
        'CYX','CYSF', 'CYSP', 'LYSR',
        'MPTYR','DPTYR','MPSER','DPSER','MPTHR','DPTHR']

    for r in aa:
        if(r == 'GLY'):
            apply_pres(FFconv,"N","GLYP", r,"N"+r)
        elif(r=='PRO'):
            apply_pres(FFconv,"N","PROP", r,"N"+r)
        else:
            apply_pres(FFconv,"N","NTER", r,"N"+r)
        apply_pres(FFconv,"C","CTER", "N"+r, "NC"+r)

    for tmp in ["CTER","CT2","CNEU"]:
        for r in aa:
            apply_pres(FFconv,"C",tmp, r,r+"_"+tmp)

    for r in aa:
        if(r in ["PRO","GLY"]): continue
        apply_pres(FFconv,"C","CT1", r,r+"_CT1")
        apply_pres(FFconv,"N","NNEU", r,"NNEU_"+r)
    return


def patch_nucleic_acids(FFconv):
    pyrimidines=["CYT","URA","THY"]
    purines=["ADE","GUA"]

    na=[]
    for r in pyrimidines:
        apply_pres(FFconv,"","DEO1", r,"D"+r) # oldstyle
        apply_pres(FFconv,"","DEOX", r,"D"+r) # newstyle
        na.extend([r,"D"+r])

    for r in purines:
        apply_pres(FFconv,"","DEO2", r,"D"+r) # oldstyle
        apply_pres(FFconv,"","DEOX", r,"D"+r) # newstyle
        na.extend([r,"D"+r])

    fiveprime =["5TER","5MET","5PHO","5POM","5DP"]
    threeprime=["3TER","3PO3","3PHO","3POM"]
    for r in na:
        for p5 in fiveprime:
            name=p5+"_"+r
            apply_pres(FFconv,"5",p5, r,name)
            for p3 in threeprime:
                apply_pres(FFconv,"3",p3, name, name+"_"+p3)
        for p3 in threeprime:
            apply_pres(FFconv,"3",p3, r,r+"_"+p3)

def patch_lipids(FFconv):

    def build_lipid(name,res1,res2,res3,patch12,patch13):
        tsys=FFconv.viparrff.typer.findTemplate(name)
        if len(tsys)>0: return
        tmpres=mergeResiduesViaPatch(FFconv, res1, res2, patch12, ("","_A"))

        if(tmpres is not None):
            tmpres.resname=name+"_TMP"
            matches,err=FFconv.viparrff.typer.matchFragment(tmpres,tmpres.system.atoms)
            assert(len(matches)==0)
            FFconv.viparrff.typer.addTemplate(tmpres)
            finalres=mergeResiduesViaPatch(FFconv, name+"_TMP", res3, patch13, ("","_B"))
            FFconv.viparrff.typer.delTemplate(tmpres)
            if(finalres is not None):
                finalres.resname=name
                ffconverter.addTemplateData(FFconv.viparrff.typer, finalres, True)

    # build DPPC from PCGL,PALM,PALM,EST1,EST2
    build_lipid("DPPC", "PCGL", "PALM","PALM","EST1","EST2")
    # build DOPC from PCGL,OLEO,OLEO,EST1,EST2
    build_lipid("DOPC", "PCGL", "OLEO","OLEO","EST1","EST2")
    # build DSPC from PCGL,STEA,STEA,EST1,EST2
    build_lipid("DSPC", "PCGL", "STEA","STEA","EST1","EST2")
    # build SDPC from PCGL,STEA,DHA,EST1,EST2
    build_lipid("SDPC", "PCGL", "STEA","DHA","EST1","EST2")



def parse_decl(FFconv,lines):
    for l in lines:
        (words,comment)=clean_line(lines[0])
        if(len(words)==0): continue
        if(len(words)!=2):
            print("Invalid DECL statement (2): %s"%(str(words)))
            continue
        decl=words[1]
        if(decl in FFconv.bsplit): continue
        i=len(FFconv.bsplit)
        FFconv.bsplit[decl]="$%d"%(i+1)

def loadparam_and_template(filename, paramonly, FFconv):
    ffsections = {
        # param sections
        ("BOND", 0) : parse_stretch,
        ("ANGL", 0) : parse_angle,
        ("THET", 0) : parse_angle,
        ("DIHE", 0) : parse_proper,
        ("PHI" , 0) : parse_proper,
        ("IMPR", 0) : parse_improper,
        ("IMPH", 0) : parse_improper,
        ("CMAP", 0) : parse_cmap,
        ("NONB", 1) : parse_lj,
        ("NBON", 1) : parse_lj,
        ("NONB", 0) : parse_lj,
        ("NBON", 0) : parse_lj,
        ("HBON", 1) : parse_nop,
        ("NBFI", 0) : parse_NBFIX,
        # template sections
        ("DECL", 1) : parse_decl,
        ("MASS", 1) : parse_mass,
        ("RESI", 1) : parse_residue,
        ("PRES", 1) : parse_residue,
        ("DEFA", 1) : parse_nop,
        ("AUTO", 1) : parse_nop,
        ("END",  0) : parse_nop
    }

    file = open(filename,'r')
    pdata=file.readlines()
    file.close()

    sections=[]
    sstop=[]
    for idx,l in enumerate(pdata):
        if len(l)==0: continue
        # strip comments and equals
        l=l.replace('=',' ').replace('!',' ! ')
        pdata[idx]=l
        words,comment=clean_line(l)
        if(len(words)==0): continue
        name=words[0][0:4].upper()
        if (len(words)>1):
            key=(name,1)
        else:
            key=(name,0)
        if key in ffsections:
            sstop.append(idx)
            while(pdata[idx].split()[-1] == '-'): idx+=1
            sections.append((key,idx))
    if(len(sstop)>0):
        sstop.pop(0)
        sstop.append(len(pdata))

    for (sn,sb),se in zip(sections,sstop):
        try:
            if paramonly:
                if (sn[0]!="RESI" and sn[0]!="PRES"):
                    ffsections[sn](FFconv,pdata[sb:se])
            else:
                if (sn[0]=="RESI" or sn[0]=="PRES"):
                    ffsections[sn](FFconv,pdata[sb:se])
        except:
            print("Failed parsing file,section,line: ",filename,sn,sb)
            raise


def setCharmmForcefieldFiles(baseMerge,parameters,templates):

    if(baseMerge!=""):
        baseMerge=os.path.abspath(baseMerge)
        if( not os.path.exists(baseMerge)):
            raise UserWarning("Couldnt find supplied base ff dir to allow merging: "+baseMerge)

    params,templ = list(),list()
    for p in parameters:
        pfull=os.path.abspath(p)
        if(not os.path.exists(pfull)): raise UserWarning("Couldnt find supplied parameter file: "+pfull)
        if(pfull not in params): params.append(pfull)
    for fftype,t in templates:
        tfull=os.path.abspath(t)
        if(not os.path.exists(tfull)): raise UserWarning("Couldnt find supplied template file: "+tfull)
        #type=getTemplateFileClassFromFileName(t,fftype)
        data=(t,"",fftype)
        if(data not in templ): templ.append(data)

    cont=Namespace()
    cont.baseMerge=baseMerge
    cont.params=params
    cont.templ=templ
    return cont

def convertCharmmForcefield(ffCharmm, nbfixID, merge, skipres, onlyThese, remove_cmap, cleanTemplates, patch_custom=None):
    cont=Namespace()
    cont.viparrff=ffconverter.initializeSimpleForcefield([0,0,1.0],[0,0,1.0],"charmm")

    def fix_atypes(_atypes):
        # Better atom type names
        conv_atype={
            "ZN"  : "Zn2+",
            "MG"  : "Mg2+",
            "CLA" : "Cl-",
            "CAL" : "Ca2+",
            "BAR" : "Ba2+",
            "CAD" : "Cd2+",
            "LIT" : "Li+",
            "SOD" : "Na+",
            "POT" : "K+",
            "RUB" : "Rb+",
            "CES" : "Cs+",
        }
        if(isinstance(_atypes,list)):
            atypes=list(map(str.upper,_atypes))
            atypes=[conv_atype.get(at,at) for at in atypes]
        else:
            atypes=_atypes.upper()
            atypes=conv_atype.get(atypes,atypes)
        return atypes
    cont.fix_atypes=fix_atypes

    #conversion factors for charmm ff units into viparr ff units

    cont.conv_lj       = (2.0 / 2.0**(1./6.), -1 ) # sig,eps
    cont.conv_nbfix_lj = (1.0 / 2.0**(1./6.), -1 ) # sig,eps
    cont.atype2atomicnumber={}
    cont.bsplit={"+N":"$1","-C":"$2"}
    cont.patches={}
    cont.skipres = [r for r in skipres if r not in onlyThese]
    cont.onlyThese=onlyThese
    cont.remove_cmap=remove_cmap
    cont.nbfixID=nbfixID

    if(merge != ''):
        preload_merge(merge,cont)

    for pload in [True, False]:
        if(pload):
            for f in ffCharmm.params:
                loadparam_and_template(f,pload,cont)

        currenttfile=''
        for file,pref,type in ffCharmm.templ:
            loadparam_and_template(file,pload,cont)
        if(not pload):
            patch_amino_acids(cont)
            patch_nucleic_acids(cont)
            patch_lipids(cont)
            if(patch_custom is not None):
                patch_custom(cont)

    typer=cont.viparrff.typer
    if(not len(onlyThese)):
        for name in skipres:
            for tmpl in typer.findTemplate(name):
                typer.delTemplate(tmpl)

    if(cleanTemplates):
        deletable=[]
        for tmpl in typer.templates:
            for atom in tmpl.system.atoms:
                if(atom.atomic_number<1): continue
                if(tmpl.btype(atom) not in cont.atype2atomicnumber or tmpl.nbtype(atom) not in cont.atype2atomicnumber):
                    print(tmpl.btype(atom),tmpl.nbtype(atom))
                    deletable.append(tmpl)
                    break
        if(len(deletable)):
            print("Found templates with missing typenames... Deleting: ",[t.name for t in deletable])
            for tmpl in deletable:
                typer.delTemplate(tmpl)

    if(remove_cmap):
        print("User requested cmap removal")
        cont.viparrff.rules.plugins=[p for p in cont.viparrff.rules.plugins if p not in ["torsiontorsion_cmap","cmap"]]
        for template in cont.viparrff.typer.templates:
            for cmapAtoms in template.cmaps:
                template.removeCmap(cmapAtoms)

    return cont.viparrff

######################################################################

def main():
    import argparse

    class TemplateGroupAppend(argparse.Action):
        def __init__(self, option_strings, dest, **kwargs):
            super(TemplateGroupAppend, self).__init__(option_strings, dest, **kwargs)

        def __call__(self, parser, namespace, values, option_string=None):
            list_object = getattr(namespace, self.dest, None) or list()
            ttype=getattr(namespace, 'tgroup', '')
            print("'%s' '%s'"%(ttype,values))
            list_object.append((ttype,values))
            setattr(namespace, self.dest, list_object)

    desc = '''
    ff_charmm_to_viparr is a force field conversion program.
    -t and -p are used to specify force field topology and parameter files
    respectivly; the order of topology and parameter files are important:
    If there are conflicts, earlier topologies/parameters take precedence over later.
    -n can be used to group templates into particular template files
    '''

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('ffdir',type=str,action='store',help="output forcefield directory")
    parser.add_argument('-n', type=str, dest='tgroup', default='',
                        help='''Group subsequent templates within this template file.
                        Multiple template groups can be specified, each with its own -n''')
    parser.add_argument('-t', type=str, action=TemplateGroupAppend, dest='templates',
                        help='''template file to parse; several can be listed,
                        each preceded by its own -t''')
    parser.add_argument('-p', action='append', dest='parameters', default=[],
                        help='''parameter file to parse; several can be listed,
                        each preceded by its own -p''')
    parser.add_argument('-r', '--remove-cmap',action='store_true', dest='remove_cmap', default=False,
                        help='''force removal of cmap terms''')
    parser.add_argument('--nbfixID',action='store',default="charmm",help="nbfix identifier for forcefield")
    parser.add_argument('-s', '--skip-res',action='append', dest='skipres', default=["TIP3","TP3M", "DUM"],
                        help='''skip residues with these names''')
    parser.add_argument('--only', type=str, action='append',dest='only', default=[],
                        help='''only keep these templates''')
    parser.add_argument('--merge',type=str, dest='merge', default='',
                        help='''base viparr directory for merging''')
    parser.add_argument('--cleanTemplates', action='store_true', default=False,
                        help='''remove templates with missing atomtypes from templates file''')

    print("CMD: "+" ".join(sys.argv))
    args = parser.parse_args()
    if(len(args.parameters)==0 and len(args.templates)==0):
        parser.error('must provide -p or -t options')

    ffCharmm = setCharmmForcefieldFiles(args.merge,args.parameters,args.templates)
    ffViparr = convertCharmmForcefield(ffCharmm, args.nbfixID,args.merge, args.skipres, args.only, args.remove_cmap, args.cleanTemplates)

    dirname=os.path.abspath(args.ffdir)
    if os.path.isdir(dirname):
        import shutil
        shutil.rmtree(dirname)

    viparr.ExportForcefield(ffViparr,dirname)

    # if(args.merge): os.unlink(os.path.join(dirname,"rules"))



