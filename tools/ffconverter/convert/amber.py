from __future__ import print_function

import os
import sys
from argparse import Namespace
from collections import OrderedDict
import msys, viparr

from viparr import ffconverter

def aname2elem(a):
    '''Given an atom name, return the element name.  This is not fool-proof.
    '''
    a = a.lstrip('1234') # remove any leading digits 1 to 4
    # assume first char and subsequent lower-case chars form element name
    name = a[0]
    for c in a[1:]:
        if c.islower():
            name = name + c
        else:
            break
    return name


import re
# try to be as forgiving as possible in parsing the files
_reSubs=dict(
    pathname=r'''(".+?"|'.+?'|[\S]+)''',                         # quoted string or non-whitespace
    fnum=r'''([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)''', # %e,%E, %f, %g compatable
    atomtype=r'''([^\s'"]+-?)''',                                 # anything but quotes,whitespace and -
    comment=r'''((?:\s+.*)|(?:[^\d\.eE+-].*))?'''               # can be space_anything or nospace_notCharsInFnum_anything
)

def makeRe(ntype,nparams):
    s='^\s*'
    if(ntype>0):
        s+="\s*-\s*".join([r'%(atomtype)s']*ntype)
        s+=r'\s+%(fnum)s'*nparams
    elif(nparams>0):
        s+="\s+".join([r'%(fnum)s']*nparams)
    s+=r'%(comment)s$'
    s=s%_reSubs
    return re.compile(s%_reSubs)

_re_bonds = makeRe(2,2)
_re_angles = makeRe(3,2)
_re_propers = [makeRe(4,4), makeRe(0,4)]
_re_impropers = [makeRe(4,4), makeRe(4,3)]
_sceere = re.compile(r'SCEE=\s*%(fnum)s' % _reSubs,re.IGNORECASE)
_scnbre = re.compile(r'SCNB=\s*%(fnum)s' % _reSubs,re.IGNORECASE)

verbose=False

def error(l):
    print('bad input at line: ', l)
    raise IOError

def sec_parse_from_re(lines, regex_list_or_single):

    if(isinstance(regex_list_or_single, (list,tuple))):
        attempts=regex_list_or_single
    else:
        attempts=[regex_list_or_single]

    parsed=[]
    for l in lines:
        lclean=l.strip()
        if(len(lclean)==0): break
        found=[]
        for regex in attempts:
            match=regex.match(lclean)
            if(match): found.append(match)
        if(len(found)==0):
            print("Failed to match regex to line: ", lclean)
            print("  Regex(s) ="," ".join([regex.pattern for regex in attempts]))
            assert(False)
        parsed.append(["" if m is None else m for m in found[0].groups()])
    return parsed

def parse_asym(file,FFconv,modify):
    pat=re.compile(makeRe(1,2).pattern+'|'+makeRe(1,1).pattern)
    for raw in sec_parse_from_re(file,pat):
        data=[d for d in raw if d != ""]
        assert(len(data) in (2,3,4))
        type=FFconv.fix_atypes(data[0])
        FFconv.seenAtomSyms.add(type)
        params={
            "type": type,
            "amu": float(data[1]),
            "memo": data[-1].strip() if len(data)>2 else ""
        }
        ffconverter.addOrUpdateParameterData(FFconv.viparrff,"mass", params, modify, False)

def parse_hydrophilic(file):
    while (1):
        l=next(file)
        nhydro=len(l.split())
        if nhydro < 20:
            break

def parse_stretch(file,FFconv,modify):
    for (at0,at1,fc,r0,memo) in sec_parse_from_re(file,_re_bonds):
        params={
            "type": " ".join(map(FFconv.fix_atypes,[at0,at1])),
            "fc"  : float(fc),
            "r0"  : float(r0),
            "memo": memo.strip()
        }
        ffconverter.addOrUpdateParameterData(FFconv.viparrff,"stretch_harm", params, modify, False)

def parse_angle(file,FFconv,modify):
    for (at0,at1,at2,fc,theta0,memo) in sec_parse_from_re(file,_re_angles):
        params={
            "type"  : " ".join(map(FFconv.fix_atypes,[at0,at1,at2])),
            "fc"    : float(fc),
            "theta0": float(theta0),
            "memo"  : memo.strip()
        }
        ffconverter.addOrUpdateParameterData(FFconv.viparrff,"angle_harm", params, modify, False)

def process_torsion(FFconv,isFrcmod,ttype,data):
    if(ttype=="proper"):
        wild=FFconv.wildprop
        nonwild=FFconv.nonwildprop
    elif(ttype=="improper"):
        wild=FFconv.wildimprop
        nonwild=FFconv.nonwildimprop
    else:
        assert(False)
    assert isinstance(wild,dict)
    assert isinstance(nonwild,dict)
    paramtype="frcmod" if(isFrcmod) else "parm"

    lastType=()
    lastCont=False
    for (atypes,fc,phase,pn,memo) in data:
        if(lastCont and lastType!=atypes):
            print('Suspect multiline %s torsion. Previous: "%s"  Current "%s"'%(ttype,lastType,type))
        lastType=atypes
        lastCont=pn<0
        pn = abs(pn)

        # We store the dihedral data temporarily so we can handle merges and overrides properly.
        mydict=wild[paramtype] if("*" in atypes) else nonwild[paramtype]
        assert isinstance(mydict,OrderedDict)

        key=tuple(atypes)
        if(key not in mydict): mydict[key]={}
        if(pn not in mydict[key]): mydict[key][pn]=[]
        mydict[key][pn].append((phase,fc,memo.strip()))

def parse_proper(file,FFconv,isFrcmod):
    data=[]
    last = None
    for r in sec_parse_from_re(file,_re_propers):
        if len(r)==9:
            atypes = list(map(FFconv.fix_atypes,r[:4]))
            idivf,pk,phase,pn = list(map(float,r[4:8]))
            comment = r[8]
            last = atypes
        elif len(r)==5:
            atypes = last
            idivf,pk,phase,pn = list(map(float,r[:4]))
            comment = r[4]
        else:
            raise RuntimeError("Unsupported r")
        print(atypes,idivf,pk,phase,pn)

        scee=_sceere.search(comment)
        if scee is None:
            scee=1.2
        else:
            comment=comment.replace(scee.group(0),'')
            scee=float(scee.group(1))

        scnb=_scnbre.search(comment)
        if scnb is None:
            scnb=2.0
        else:
            comment=comment.replace(scnb.group(0),'')
            scnb=float(scnb.group(1))
        if(scee!=1.2 or scnb!=2.0):
            print("Unsupported SCEE or SCNB parameters in proper torsion spec: SCEE=%f SCNB=%f"%(scee,scnb))
        assert(idivf.is_integer() and pn.is_integer())
        pn=int(pn)
        if(abs(pn)>6 or abs(pn)<1):
            raise RuntimeError("Proper torsion phase is unphysical in file %s: %d"%(file,pn))
        data.append([atypes,pk/idivf,phase,pn,comment.strip()])
    process_torsion(FFconv,isFrcmod,"proper",data)


def parse_improper(file,FFconv,isFrcmod):
    data=[]
    for r in sec_parse_from_re(file,_re_impropers):
        new=list(map(FFconv.fix_atypes,r[:4]))
        rest=r[4:]
        if(len(rest)==5):
            idivf,fc,phase,pn=list(map(float,rest[:-1]))
        elif(len(rest)==4):
            idivf=1.0
            fc,phase,pn=list(map(float,rest[:-1]))
        assert(idivf.is_integer() and pn.is_integer())
        pn=int(pn)
        if(abs(pn)>6 or abs(pn)<1):
            raise "Improper torsion phase is unphysical in file %s: %d"%(file,pn)
        data.append([new,fc,phase,pn,rest[-1]])
    process_torsion(FFconv,isFrcmod,"improper",data)

def parse_10_12(file,modify):
    for l in file:
        words=l.split()
        if len(words)==0:
            return
        a1=words[0]
        a2=words[1]
        A =float(words[2])
        B =float(words[3])
        comment=words[4:]

        if A!=0.0 or B!=0.0:
            print("10-12 potential conversion not supported at this time")
            error(l)

    return

def parse_equiv_6_12(file,FFconv):
    for l in file:
        words=l.split()
        if len(words)==0:
            return

        key=FFconv.fix_atypes(words)
        org =key[0]
        equiv=key[1:]
        for t in equiv:
            FFconv.nbequiv[t]=org

    return

def parse_6_12(file,FFconv,modify):
    if(not modify):
        l=next(file)
        words=l.split()
        if len(words) != 2:
            error(l)
        type=words[1]

        if type == 'RE':
            if FFconv.vdwtype is not None: assert(FFconv.vdwtype=="se")
            FFconv.vdwtype = 'se'
        elif type == 'AC':
            if FFconv.vdwtype is not None: assert(FFconv.vdwtype=="cc")
            FFconv.vdwtype = 'cc'
        else:
            print("Unsupported Type")
            error(l)

    if(FFconv.vdwtype is None): FFconv.vdwtype = 'se'
    if (FFconv.vdwtype=='se'):
        conv = FFconv.conv_lj_se
    else:
        raise UserWarning("Unsupported Type: "+FFconv.vdwtype)


    for l in file:
        words=l.split()
        if(len(words) ==0):
            if(not modify):
                continue
            else:
                return
        if( words[0].upper() == 'END'):
            return
        key = FFconv.fix_atypes(words[0])
        params={
          "type": key,
          "sigma": float(words[1]),
          "epsilon": float(words[2]),
          "memo":" ".join(words[3:]),
          "nbfix_identifier":""
        }
        for idx, pname in enumerate(["sigma","epsilon"]):
            params[pname]*=conv[idx]

        ffconverter.addOrUpdateParameterData(FFconv.viparrff,"vdw1", params, modify, False)

    return

def parse_cmap(file, FFconv, modify):
    allcmaps = []
    cmap = {}
    mode = ""
    nmode = None

    for l in file:
        lclean=l.strip()
        if(len(lclean)==0): break
        words = lclean.split()
        if words[0] == "%FLAG":
            if mode not in ["", "CMAP_TITLE"]:
                assert len(cmap[mode]) == nmode
            elif mode == "CMAP_TITLE":
                cmap[mode] = " ".join(cmap[mode])

            mode = ""
            nmode = None
            key = words[1]

            if key == "CMAP_TITLE":
                if len(cmap):
                    allcmaps.append(cmap)
                    cmap = {}
            assert key not in cmap
            # Its unclear from seeing two different files done differently if this is the
            # index of the CMAP table or how many cmap tables to expect... funtimes
            if key == "CMAP_COUNT":
                pass
            elif key == "CMAP_TITLE":
                mode = key
                nmode = None
                cmap[mode] = []
            elif key == "CMAP_RESLIST":
                mode = key
                nmode = int(words[2])
                cmap[mode] = []
            elif key == "CMAP_RESOLUTION":
                cmap[key] = int(words[2])
            elif key == "CMAP_ATMLIST":
                cmap[key] = words[2].split('-')
            elif key == "CMAP_RESIDX":
                cmap[key] = list(map(int, words[2:]))
            elif words[1] == "CMAP_PARAMETER":
                mode = "CMAP_PARAMETER"
                nmode = cmap["CMAP_RESOLUTION"]**2
                cmap[mode] = []
            else:
                raise RuntimeError(f"Unknown CMAP flag: '{key}'")
        elif(words[0] == "%COMMENT"):
            cmap["COMMENT"] = cmap.get("COMMENT", "") + " " + " ".join(words[2:])
        elif mode in ["CMAP_RESLIST", "CMAP_PARAMETER", "CMAP_TITLE"]:
            cmap[mode].extend(words)
        else:
            raise RuntimeError("Unhandled file structure")

    if len(cmap):
        allcmaps.append(cmap)

    seen = set()
    for cmap in allcmaps:
        cmapTable = msys.CreateParamTable()
        cmapTable.addProp('phi', float)
        cmapTable.addProp('psi', float)
        cmapTable.addProp('energy', float)
        FFconv.viparrff.addCmapTable(cmapTable)

        spacing = 360.0/cmap["CMAP_RESOLUTION"]
        evals = list(map(float, cmap["CMAP_PARAMETER"]))
        for iphi in range(cmap["CMAP_RESOLUTION"]):
            for ipsi in range(cmap["CMAP_RESOLUTION"]):
                param = cmapTable.addParam()
                param["phi"] = -180.0 + iphi*spacing
                param["psi"] = -180.0 + ipsi*spacing
                param["energy"] = evals[cmap["CMAP_RESOLUTION"] * iphi + ipsi]

    if len(FFconv.cmaps) and len(allcmaps):
        raise RuntimeError("I cant overwrite cmap tables yet")
    FFconv.cmaps = allcmaps

# FIXME: This is hacked together to support the amber ff19SB forcefield
# with our current viparr code.
# I have no illusions that it will work for general cmap terms in amber
# which would require modifications of viparr
def apply_cmap_to_templates(FFconv):
    # we are going to go through, and rename the bonded type used for the
    # alpha carbon of each cmap term, and then duplicate the bonded parameters as
    # necessary.
    # Also make sure all residues use the same atomtype names
    if not FFconv.cmaps: return
    unique_types = set()
    bbnames =  ["N", "CA", "C"]
    typer = FFconv.viparrff.typer
    rewrite = []
    for icmap, cmap in enumerate(FFconv.cmaps):
        assert "CMAP_RESIDX" not in cmap
        assert "CMAP_ATMLIST" not in cmap
        title = cmap["CMAP_TITLE"]

        for resname in cmap["CMAP_RESLIST"]:
            tmpls = typer.findTemplate(resname)
            assert len(tmpls)
            for tmpl in tmpls:
                found = []
                for aname in bbnames:
                    for a in tmpl.system.atoms:
                        if a.name == aname:
                            found.append((a, tmpl.btype(a)))
                            break
                assert len(found) == len(bbnames)
                unique_types.add(tuple([v[1] for v in found]))
                assert len(unique_types)==1
                ca, oldtype = found[1]
                newtype = oldtype+"_%02d"%(icmap+1)
                rewrite.append((oldtype, newtype))
                tmpl.setTypes(ca, newtype, tmpl.nbtype(ca), tmpl.pset(ca))
                found = [v[0] for v in found]
                am1 = [a for a in found[0].bonded_atoms if a.atomic_number==-1]
                ap1 = [a for a in found[-1].bonded_atoms if a.atomic_number==-1]
                assert(len(am1)==1 and len(ap1)==1)
                cmapatoms = am1 + found + found + ap1
                tmpl.addCmap(cmapatoms)
        atypes = list(list(unique_types)[0])
        atypes[1] = rewrite[-1][1]
        t1 = ["*"]+atypes
        t2 = atypes+["*"]
        params = {
                "type":" ".join(t1+t2),
                "cmapid": icmap+1,
                "memo": title
                }
        print(params)
        ffconverter.addOrUpdateParameterData(FFconv.viparrff,"torsiontorsion_cmap",params,False, False, False)
    if len(rewrite)==0: return
    rewrite = set(rewrite)
    rewriteRef = set([v[0] for v in rewrite])
    assert len(rewriteRef) == 1, f"Got rewriteRef {rewriteRef}; expected exactly 1 element"
    rewriteRef = rewriteRef.pop()
    rewrite = sorted([v[1] for v in rewrite])
    assert len(set(rewrite)) == len(rewrite)
    for tname in ['mass', 'stretch_harm', 'angle_harm', 'dihedral_trig', 'improper_trig']:
        table = FFconv.viparrff.ParamTable(tname)
        params = table.params
        FFconv.viparrff.delParams(tname, params=params)
        for p in params:
            kwds = { k:p[k] for k in p.keys() }
            FFconv.viparrff.appendParam(tname, **kwds)
            atypes = kwds["type"].split(" ")
            if rewriteRef in atypes:
                kwds["memo"] += " VIPARR CMAP SUPPORT"
                for rw in rewrite:
                    kwds["type"] = " ".join([at if at != rewriteRef else rw for at in atypes])
                    FFconv.viparrff.appendParam(tname, **kwds)

        #print(params[0].items())






# amber handles dihedral parameter overides differently than most forcefields.
# Essentially, terms only update the particular multiplicity they refer too, even
# if they occur at a later point in the file (the -pn  spec is completly bogus)
# wild card terms are applied first, followed by specific terms. Thus, if you dont include
# all the phases for the non-wild term, you will get a combination of wild/nonwild parameters
# If you want to update a torsion later (eg frcmod) you need to specify all multiplicites or again
# you can end up with part old parameters and part new... Not very pretty...
# WOW - its even worse than I initially thought. You DONT use the wild card dihedral terms
# for the specific terms if you are in a frcmod file... what a mess...
def merge_and_add_torsions(FFconv,name):
    if(name=="dihedral_trig"):
        wild=FFconv.wildprop
        nonwild=FFconv.nonwildprop
    elif(name=="improper_trig"):
        wild=FFconv.wildimprop
        nonwild=FFconv.nonwildimprop
    else:
        assert(False)

    assert isinstance(wild,dict)
    assert isinstance(wild["parm"],OrderedDict)
    assert isinstance(wild["frcmod"],OrderedDict)
    assert isinstance(nonwild,dict)
    assert isinstance(nonwild["parm"],OrderedDict)
    assert isinstance(nonwild["frcmod"],OrderedDict)

    zero=0.0

    merged=OrderedDict(wild["parm"])
    tmpmerge=OrderedDict()
    # make sure nonwild dihedrals contain all necessary terms
    # (inherited from wilds BUT only for non-frcmod parameters)
    for k in nonwild["parm"]:
        atypesFwd=("*",k[1],k[2],"*")
        atypesRev=tuple(atypesFwd[::-1])
        if(name=="dihedral_trig" and atypesFwd!=atypesRev):
            check=[atypesFwd,atypesRev]
        else:
            check=[atypesFwd]
        for atypes in check:
            found=0
            if atypes in wild["parm"]:
                found+=1
                for m in wild["parm"][atypes]:
                    if m in nonwild["parm"][k]: continue
                    nonwild["parm"][k][m]=wild["parm"][atypes][m]
            assert(found<2)
        tmpmerge[k]=nonwild["parm"][k]

    # Now append frcmod parameters
    for tmp in [wild["frcmod"],tmpmerge,nonwild["frcmod"]]:
        for keyFwd,value in tmp.items():
            keyRev=tuple(keyFwd[::-1])
            if(name=="dihedral_trig" and keyFwd!=keyRev):
                assert(not (keyFwd in merged and keyRev in merged))
                if(keyRev in merged):
                    key=keyRev
                else:
                    key=keyFwd
            else:
                key=keyFwd
            if key not in merged: merged[key]={}
            for pn,data in value.items():
                if pn not in merged[key]: merged[key][pn]=[]
                merged[key][pn].extend(data)

    for k in merged:
        tlist=[]
        for pn in merged[k]:
            # take the last match for a given multiplicity
            phase,fc,comment=merged[k][pn][-1]
            idxlist=range(len(tlist))
            data=[zero]*8

            if len(idxlist)==0:
                data[0]=zero
                data[1]=fc
                if (phase==zero):
                    data[1+pn]=fc
                elif (phase==180.0):
                    data[1+pn]=-fc
                else:
                    data[0]   =phase
                    data[1+pn]=fc
                tlist.append([data,comment])
            else:
                found=False
                for idx in idxlist:
                    data=tlist[idx][0]
                    if( (phase == zero or phase == 180.0) and data[0] == zero ):
                        found=True
                        if data[1+pn] != zero: # Should never trigger:
                            raise UserWarning("multiplicity has already been assigned a value: %s %d %f %f"%(str(atypes),pn,data[1+pn],fc))
                        data[1]+=fc
                        if (phase==zero):
                            data[1+pn]=fc
                        else:
                            data[1+pn]=-fc
                        break
                    # merge for phase != {0,180}
                    elif(phase == data[0]):
                        found=True
                        if data[1+pn] != zero: # Should never trigger
                            raise UserWarning("multiplicity has already been assigned a value: %s %d %f %f"%(str(atypes),pn,data[1+pn],fc))
                        data[1]+=fc
                        data[1+pn]=fc
                        break
                if(found):
                    tlist[idx][0]=data
                    if(comment.strip() not in tlist[idx][1] ): tlist[idx][1]=tlist[idx][1]+" "+comment.strip()
                else:
                    data=[zero]*8
                    data[0]=phase
                    data[1]+=fc
                    data[1+pn]=fc
                    tlist.append([data,comment])

        if(name=="dihedral_trig"):
            groupedList=[]
            for t in tlist:
                order=["phi0","fc0","fc1","fc2","fc3","fc4","fc5","fc6"]
                assert len(t[0])==len(order)
                params=dict(zip(order,t[0]))
                params["memo"]=t[1]
                params["type"]=" ".join(k)

                if(len(groupedList) and groupedList[-1][0]["type"]==params["type"]):
                    groupedList[-1].append(params)
                else:
                    groupedList.append([params])
            for params in groupedList:
                ffconverter.addOrUpdateParameterData(FFconv.viparrff, name, params, False, True)
        else:
            for t in tlist:
                order=["phi0","fc0","fc1","fc2","fc3","fc4","fc5","fc6"]
                assert len(t[0])==len(order)
                params=dict(zip(order,t[0]))
                params["memo"]=t[1]
                params["type"]=" ".join(k)
                ffconverter.addOrUpdateParameterData(FFconv.viparrff, name, params, False, False)


def load_frcmod_or_param(filename, FFconv):
    file = open(filename,'r')

    # Section 1 Title
    title=next(file)

    datatype=None
    for line in file:
        l=line.strip()
        if(len(l)==0): continue

        if (l.split()[0].upper()[:4] in ['MASS','BOND','ANGL','DIHE','IMPR','HBON','NONB']):
            datatype="frcmod"
        else:
            datatype="param"
        break
    file.close()

    if(datatype is None):
        print("Nothing to load in file: ",filename)
        return
    print("loading: %s  filetype: %s"%(filename,datatype))
    if(datatype=="frcmod"):
        _loadfrcmod(filename,FFconv)
    else:
        _loadparam(filename,FFconv)
    print("done loading: ",filename)

def _loadfrcmod(filename,FFconv):
    file = open(filename, 'r')

    # Section 1 Title
    l=next(file)

    for l in file:
        words=l.strip().split()
        if(len(words) ==0): continue

        if words[0][0:4] == 'MASS':
            parse_asym(file,FFconv,True)
        elif words[0][0:4] == 'BOND':
            parse_stretch(file,FFconv,True)
        elif words[0][0:4] == 'ANGL':
            parse_angle(file,FFconv,True)
        elif words[0][0:4] == 'DIHE':
            parse_proper(file,FFconv,True)
        elif words[0][0:4] == 'IMPR':
            parse_improper(file,FFconv,True)
        elif words[0][0:4] == 'HBON':
            parse_10_12(file,True)
        elif words[0][0:4] == 'NONB':
            parse_6_12(file,FFconv,True)
        elif words[0][0:4] == 'CMAP':
            parse_cmap(file, FFconv, True)
    file.close()

def _loadparam(filename,FFconv):
    file = open(filename,'r')

    # Section 1 Title
    l=next(file)

    # Section 2 Atom Symbols and masses
    parse_asym(file,FFconv,False)

    # Section 3 Hydrophobic symbols
    parse_hydrophilic(file)

    # Section 4 bond parameters
    parse_stretch(file,FFconv,False)

    # Section 5 angle
    parse_angle(file,FFconv,False)

    # Section 6 dihedral
    parse_proper(file,FFconv,False)

    # Section 7 impropers
    parse_improper(file,FFconv,False)

    # Section 8 10-12 vdw parameters
    parse_10_12(file,False)

    # Section 9 equivalence for 6-12 parameters
    parse_equiv_6_12(file,FFconv)

    # Section 10 6-12 parameters
    parse_6_12(file,FFconv,False)

    file.close()
    return


# This is expecting a *.in file. It should be
# extended to support mol2 at some point
# (amber seems to be moving that direction)
def loadtempl(filename,pref,fftype,FFconv):
    print("LOADING TMPL", filename)
    file = open(filename,'r')

    next(file) # card1
    next(file) # card2

    while(1):
        words=next(file).split()
        #print words
        if len(words)==0:
            break
        if words[0] == 'STOP':
            break

        #card 3
        title=' '.join(words)

        #card 4
        namef=next(file)

        #card 5
        tname=next(file).split()[0]
        nameres=pref+tname

        #card 6
        ifixc=next(file).split()[0]

        #card 7
        next(file)

        #card 8
        tsys = viparr.TemplatedSystem()
        res = tsys.system.addResidue()
        res.name=nameres
        amap={}
        alist=[]
        mainchain=[]
        while(1):
            line=next(file).strip()
            if len(line)==0: break
            words=line.split()

            idx=int(words[0]) - 4
            if idx <0: continue
            aname=words[1]
            atype=FFconv.fix_atypes(words[2])
            atypenb=atype
            branch=words[3]
            ibond = int(words[4]) - 4
            acharge=float(words[10]) if len(words)==11 else 0.0

            a=res.addAtom()
            a.name=aname
            a.atomic_number=msys.ElementForAbbreviation(aname2elem(aname))
            a.charge=acharge
            tsys.setTypes(a,atype,atypenb)
            amap[aname]=a
            alist.append(a)

            if(branch=='M'): mainchain.append(a)

            if( ibond >=0 ):
                if(ibond==idx):
                    if(branch!='M'):
                        print("Bad bond: %s %s %s has a self connection. Reconnecting to %s: "%(os.path.basename(filename), res, a.name,mainchain[-1].name))
                        a.addBond(mainchain[-1])
                    else:
                        assert(False)
                else:
                    bonded=alist[ibond]
                    a.addBond(bonded)

        if(tsys.system.nbonds>2):
            if(fftype == "amino_acids"):
                print((fftype, pref, nameres))
                if( nameres in ['NME','NHE'] or (pref == "" and nameres not in ['ACE']) or pref == "C"):
                    a=res.addAtom()
                    a.name="$1"
                    a.atomic_number=-1
                    tsys.setTypes(a,"C","C")
                    amap[a.name]=a
                    a.addBond(mainchain[0])
                if(nameres=='ACE' or (pref == "" and nameres not in ['NME','NHE']) or pref == "N"):
                    a=res.addAtom()
                    a.name="$2"
                    a.atomic_number=-1
                    tsys.setTypes(a,"N","N")
                    amap[a.name]=a
                    a.addBond(mainchain[-1])
                if(pref not in ["C", "N", ""]):
                    print("What type of prefix is '%s'" %(pref))
                if('CYX' in tname):
                    a=res.addAtom()
                    a.name="$3"
                    a.atomic_number=-1
                    tsys.setTypes(a,"S","S")
                    amap[a.name]=a
                    a.addBond(amap["SG"])
            elif(fftype == "carb"):
                carbon_map = { "0": [], "1": [1], "2": [2], "3": [3], "4": [4], "6": [6]}
                carbon_map.update({"Z":[2,3], "Y":[2,4], "X":[2,6], "W":[3,4], "V":[3,6], "U":[4,6]})
                carbon_map.update({"T":[2,3,4], "S":[2,3,6], "R":[2,4,6], "Q":[3,4,6], "P":[2,3,4,6]})
                for rname in list(carbon_map.keys()):
                    cvals = [f"C{i}" for i in carbon_map[rname]]
                    carbon_map[rname] = cvals
                wanted = ["C1", "C2", "C3", "C4", "C5", "C6"]
                carbons = {a.name:a for a in tsys.system.atoms if a.name in wanted}
                cap = carbon_map.get(nameres[0],[])
                for icap, cname in cap:
                    atom = carbons[cname]
                    anew = res.addAtom()
                    anew.name="$0"+cname[-1]
                    anew.atomic_number=-1
                    amap[anew.name] = anew
                    anew.addBond(atom)

            elif(fftype == 'nucleic_acids'):
                raise RuntimeError("Please Fix This Section")
                if(len(tname)==3):
                    if(tname[-1] == '3'):
                        a=res.addAtom()
                        a.name="$1"
                        a.atomic_number=-1
 # FIXME                tsys.setTypes(a,"C","C")
                        amap[a.name]=a
                        a.addBond(mainchain[0])
                    elif(tname[-1]== '5'):
                        a=res.addAtom()
                        a.name="$2"
                        a.atomic_number=-1
 # FIXME                tsys.setTypes(a,"N","N")
                        amap[a.name]=a
                        a.addBond(mainchain[-1])
                    elif(tname[-1]=='N'):
                        pass
                    else:
                        print("what type of resname is '%s'" % (tname))
                else:
                    a=res.addAtom()
                    a.name="$1"
                    a.atomic_number=-1
 # FIXME            tsys.setTypes(a,"C","C")
                    amap[a.name]=a
                    a.addBond(mainchain[0])
                    a=res.addAtom()
                    a.name="$2"
                    a.atomic_number=-1
# FIXME                tsys.setTypes(a,"N","N")
                    amap[a.name]=a
                    a.addBond(mainchain[-1])
            elif(fftype in ["ions","ligand"]):
                pass
            else:
                print("what type of fftype is '%s'" %(fftype))

 #       if(tname== 'OHE' and fftype == 'nucleic_acids'):
 #           res.add_bondedterm("bonds",["H","O"])

        if ifixc.startswith('CHA'):
            next(file)

        type=0
        xtraq=list()
        xtrabond=list()
        xtraimpr=list()
        while(1):
            words=next(file).split()
            if len(words)==0:
                continue
            if words[0] == 'DONE':
                break
            elif words[0] == 'CHARGE':
                i=0
                while(1):
                    words=next(file).split()
                    if len(words)==0:
                        break
                    for q in words:
                        atom=tsys.system.atom(i)
                        atom.charge=float(q)
                        i=i+1
            elif words[0] == 'LOOP':
                while(1):
                    words=next(file).split()
                    if len(words)==0:
                        break
                    a0,a1=list(map(amap.get,words))
                    a0.addBond(a1)
            elif words[0] == 'IMPROPER':
                # Dont process impropers anymore, just use automated method to generate
                while(1):
                    if(len(next(file).strip())==0): break
            else:
                print("Unknown section in card 9")
                error(str(words))

        # Add impropers based on amber craziness....
        OMG_do_crazy_amber_impropers(FFconv,tsys,res)

        if(tsys.name not in FFconv.skipres and (len(FFconv.onlyThese)==0 or tsys.name in FFconv.onlyThese)):
            ffconverter.addTemplateData(FFconv.viparrff.typer, tsys, True)

    file.close()

def OMG_do_crazy_amber_impropers(FFconv,tsys,res):
    import itertools
    for atom in res.atoms:
        bonded=atom.bonded_atoms
        if(len(bonded)!=3 or 0 in [a.atomic_number for a in bonded]):
            continue
        avail=[(tsys.btype(a),a.id,a.name,a) for a in bonded]
        for nkeep in range(3,0,-1):
            allFound=[]
            for i,attempt in enumerate(itertools.combinations(list(range(3)),nkeep)):
                used=sorted([avail[i] for i in attempt])
                unused=sorted([('*',avail[i][1],avail[i][2],avail[i][3]) for i in set(range(3))-set(attempt)])
                combined=unused+used
                best_guess=[combined[0],combined[1],(tsys.btype(atom),atom.id,atom.name,atom),combined[2]]
                # print "Testing: ",res.name, atom.name, nkeep,i,best_guess
                found=FFconv.viparrff.findParams("improper_trig",type=" ".join([b[0] for b in best_guess]))
                assert(len(found)<2)
                if(len(found)==1):
                    allFound.append((found[0],best_guess))
            if(len(allFound)>0):
                psets=set([p[0] for p in allFound])
                if(len(psets)>1):
                    print("Multiple matches for improper: ", res.name, atom.name, allFound)
                    assert(False)
                else:
                    best_guess=allFound[0][1]
                    # print "Adding improper: %s %s  name=%s  types=%s"%(res.name, atom.name, str([b[2] for b in best_guess]),str([b[0] for b in best_guess]))
                    tsys.addImproper([b[3] for b in best_guess])
                    break

        else:
            print("Missing Improper: ",res.name, atom.name,avail)

def loadOff(filename,pref,fftype,FFconv):
    print("LOADING OFF: ",filename,fftype)
    import parmed
    templates=parmed.load_file(filename)
    for k,t in templates.items():
        if(not hasattr(t,"atoms")): continue
        if(k in FFconv.skipres or (len(FFconv.onlyThese) and k not in FFconv.onlyThese)): continue

        #t.to_structure()
        tsys = viparr.TemplatedSystem()
        res = tsys.system.addResidue()
        res.name=k
        amap={}
        for _a in t.atoms:
            assert(_a.name not in amap)
            a=res.addAtom()
            a.name=_a.name
            a.atomic_number=_a.atomic_number
            a.charge=_a.charge
            atype=FFconv.fix_atypes(_a.type)
            atypenb=FFconv.nbequiv.get(atype,atype)
            tsys.setTypes(a,atype,atypenb)
            amap[a.name]=a

        for b in t.bonds:
            a0=amap[b.atom1.name]
            a1=amap[b.atom2.name]
            a0.addBond(a1)

        if t.head is not None:
            a=res.addAtom()
            a.name="$1"
            a.atomic_number=-1
            if(fftype=="amino_acids"):
                tsys.setTypes(a,"C","C")
            amap[a.name]=a
            a.addBond(amap[t.head.name])
        if t.tail is not None:
            a=res.addAtom()
            a.name="$2"
            a.atomic_number=-1
            if(fftype=="amino_acids"):
                tsys.setTypes(a,"N","N")
            amap[a.name]=a
            a.addBond(amap[t.tail.name])
        for i,_a in enumerate(t.connections):
            a=res.addAtom()
            a.name="$%d"%(i+3)
            a.atomic_number=-1
            if(fftype=="amino_acids" and 'CYX' in res.name):
                tsys.setTypes(a,"S","S")
            amap[a.name]=a
            a.addBond(amap[_a.name])

        # Add impropers based on amber craziness....
        OMG_do_crazy_amber_impropers(FFconv,tsys,res)

        ffconverter.addTemplateData(FFconv.viparrff.typer, tsys, True)



def getTemplateFileClassFromFileName(fname,default=None):
    lcname=os.path.basename(fname).lower()
    if(default is not None):
        type=default
    elif(lcname.startswith('solv')):
        type="solvents"
    elif("nuc" in lcname):
        type="nucleic_acids"
    elif("ion" in lcname):
        type="ions"
    else:
        type="amino_acids"
    return type

def loadamberRC(filename, parmdirs, libdirs, prepdirs, skipmissing, skipsolvents):
    print("Using rc file: "+filename)

    #  regexes
    atomtypere = re.compile(r"""{\s*["']([\w\+\-]+)["']\s*["'](\w+)["']\s*"""
                           r"""["'](\w+)["']\s*}""")
    loadprepre = re.compile(r'loadamberprep\s+%(pathname)s' % _reSubs, re.I)
    loadparamsre = re.compile(r'loadamberparams\s+%(pathname)s' % _reSubs, re.I)
    loadoffre = re.compile(r'loadoff\s+%(pathname)s' % _reSubs, re.I)
    loadmol2re = re.compile(r'(\S+)\s*=\s*loadmol[23]\s*%(pathname)s' % _reSubs, re.I)

    params=[]
    templ=[]

    file = open(filename,'r')
    strings=""
    with open(filename,'r') as file:
        for l in file.readlines():
            l=l.strip()
            if(len(l)==0 or l[0]=='#'): continue
            strings+=l+'\n'

    for atype,elem,hyb in atomtypere.findall(strings):
        if(hyb=='sp2'):
            #print "SP2 Atype, elem, hyb: ",atype,elem,hyb
            pass

    for par in loadparamsre.findall(strings):
        for pdir in parmdirs:
            d=os.path.join(pdir,par)
            if(not os.path.exists(d)): continue
            print("Using parameter file: ", d)
            params.append(d)
            break
        else:
            if(skipmissing):
                print("WARNING! Couldnt find param file. Skipping: ",par)
            else:
                raise UserWarning("didnt find viable location for param file: %s\nSearch dirs:\n%s"%(par,"\n".join(parmdirs)))

    for searchDirs,matches in [(libdirs, loadoffre.findall(strings)),(prepdirs,loadprepre.findall(strings))]:
        for name in matches:
            path=None
            for ldir in searchDirs:
                d=os.path.join(ldir,name)
                if(not os.path.exists(d)): continue
                print("Using template file: ", d)
                path=d
                break
            else:
                if(skipmissing):
                    print("WARNING! Couldnt find templates... Skipping: ",name)
                    continue
                else:
                    raise UserWarning("didnt find viable location for off file: %s"%(name))

            type=getTemplateFileClassFromFileName(path)
            if(type=="solvents" and skipsolvents):
                # solvents cause problems
                print("Requested: Skipping solvent file: ",path)
                continue

            pre=''
            if(not path.endswith('.lib')):
                lcname=os.path.basename(path).lower()
                if(lcname.find('aminoct')>=0 ):
                    pre='C'
                elif (lcname.find('aminont')>=0 ):
                    pre='N'

            templ.append((d,pre,type))


    mol2s=loadmol2re.findall(strings)
    if(len(mol2s)): print("MOL2 Files: ",mol2s)

    return params,templ


def convertAmberForcefield(amberFF, skipres, onlyThese, cleanTemplates):
    cont=Namespace()
    cont.viparrff=ffconverter.initializeSimpleForcefield([0,0,0.8333], [0, 0, 0.5])

    def fix_atypes(atypes):
        # Viparr uses '*' as a wildcard. Change to non-wildcard character, and change 'X' to wild
        conv_atype={"C*" : "C&",
                    "N*" : "N&",
                    "X" : "*"
        }
        if(isinstance(atypes, list)):
            atypes=[conv_atype.get(at, at) for at in atypes]
        else:
            atypes=conv_atype.get(atypes, atypes)
        return atypes

    cont.fix_atypes=fix_atypes
    cont.cmaps = []

    cont.wildprop={"parm" : OrderedDict(), "frcmod" : OrderedDict()}
    cont.nonwildprop={"parm" : OrderedDict(), "frcmod" : OrderedDict()}
    cont.wildimprop={"parm" : OrderedDict(), "frcmod" : OrderedDict()}
    cont.nonwildimprop={"parm" : OrderedDict(), "frcmod" : OrderedDict()}

    cont.vdwtype = 'se'
    cont.nbequiv = dict()
    cont.seenAtomSyms=set()

    # conversion factors for amber ff file parameters into viparr ff file parameters
    cont.conv_lj_se    = (2.0 / 2.0**(1./6.), 1.0)  # sig,eps
    cont.conv_lj_cc    = (1.0,-1.0)              # c12,c6
    cont.skipres = skipres
    cont.onlyThese=onlyThese

    for file in amberFF.params:
        load_frcmod_or_param(file, cont)
    merge_and_add_torsions(cont, "dihedral_trig")
    merge_and_add_torsions(cont, "improper_trig")

    for file in amberFF.nbequiv:
        nstart=len(cont.nbequiv)
        parse_equiv_6_12(file, cont)
        if(len(cont.nbequiv)==nstart):
            raise RuntimeError("no entries were parsed in the extra nbequiv file... Are you sure everything is ok")

    # Add nbequiv to vdw file
    tname='vdw1'
    for new, ref in cont.nbequiv.items():
        params = cont.viparrff.findParams(tname, type=ref)
        if(len(params) == 0): raise RuntimeError("Couldnt make vdw type %s equivalent to %s. Type %s was not found"%(new,ref,ref))
        pnew = viparr.Forcefield.ParamTable(tname).addParam()
        cont.viparrff.appendParam(tname, pnew)
        for k in params[0].keys():
            pnew[k] = params[0][k]
        pnew['type'] = new
        pnew['memo'] += ' --> nbequiv to %s'%(ref)

    for t in amberFF.templ:
        file, pref, ttype = t
        s="templ: " + file + " pref: " + pref + " type:" + ttype
        print(s)
        if(file.endswith(".lib")):
            loadOff(file, pref, ttype, cont)
        else:
            loadtempl(file, pref, ttype, cont)

    apply_cmap_to_templates(cont)


    typer=cont.viparrff.typer
    # Remove water and useless Pseudo Ions
    if(not len(onlyThese)):
        for name in skipres:
            for tmpl in typer.findTemplate(name):
                typer.delTemplate(tmpl)

    if(cleanTemplates):
        print(cont.seenAtomSyms)
        deletable=[]
        for tmpl in typer.templates:
            for atom in tmpl.system.atoms:
                if(atom.atomic_number<1): continue
                if(tmpl.btype(atom) not in cont.seenAtomSyms or tmpl.nbtype(atom) not in cont.seenAtomSyms):
                    deletable.append(tmpl)
                    break
        print("Found templates with missing typenames... Deleting: ",[t.name for t in deletable])
        for tmpl in deletable:
            typer.delTemplate(tmpl)
    return cont.viparrff

def setAmberForcefieldFiles(rcpath, includes, parameters, templates, fftype, nbequiv, skipmissing, skipsolvents):
    cont=Namespace()
    if(rcpath is not None):
        ambdir=[]
        if('/cmd/' in rcpath):
            ambdir.append(rcpath.split('/cmd/')[0])
        else:
            ambdir.append(os.path.dirname(rcpath))
        for d in includes:
            if(d in ambdir): continue
            if(os.path.exists(d)): ambdir.append(d)
        print("ambdirs = ",ambdir)
        parmdir  = [ d for d in [ os.path.join(a, 'parm')       for a in ambdir ] if(os.path.exists(d)) ]
        libdir   = [ d for d in [ os.path.join(a, 'lib')        for a in ambdir ] if(os.path.exists(d)) ]
        libdir  += [ d for d in [ os.path.join(a, 'lib/oldff')  for a in ambdir ] if(os.path.exists(d)) ]
        prepdir  = [ d for d in [ os.path.join(a, 'prep')       for a in ambdir ] if(os.path.exists(d)) ]
        prepdir += [ d for d in [ os.path.join(a, 'prep/oldff') for a in ambdir ] if(os.path.exists(d)) ]

        if(len(parmdir)==0): raise UserWarning("Couldnt find any viable parm directories: searched %s"%(str(ambdir)))
        if(len(libdir)==0): raise UserWarning("Couldnt find any viable lib directories: searched %s"%(str(ambdir)))
        if(len(prepdir)==0): raise UserWarning("Couldnt find any viable prep directories: searched %s"%(str(ambdir)))

        assert(os.path.exists(rcpath))
        cont.rcfile=os.path.basename(rcpath)

        params,templ=loadamberRC(rcpath, parmdir, libdir, prepdir, skipmissing, skipsolvents)
    else:
        cont.rcfile = None
        params,templ = list(), list()

    for p in parameters:
        if(not os.path.exists(p)): raise UserWarning("Couldnt find supplied parameter file: "+p)
        params.append(p)
    for t in templates:
        if(not os.path.exists(t)): raise UserWarning("Couldnt find supplied template file: "+t)
        ttype = getTemplateFileClassFromFileName(t, fftype)
        templ.append((t, "", ttype))

    cont.params=params
    cont.templ=templ
    cont.nbequiv=nbequiv

    return cont

def main():
    import argparse

    desc = '''
    ff_amber_to_viparr is a force field conversion program.
    It has been tested against leaprc.ff{94,96,99,99SB,03} and
    may not work with other versions of the amber forcefields.
    *USER BEWARE*
  '''
    use_extra_nbequiv_files = False

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('ffdir', type=os.path.realpath, action='store', help="output forcefield directory")
    parser.add_argument('-r', type=os.path.realpath, dest='rcpath', default=None,
                         help='''path to leaprc file to convert''')
    parser.add_argument('-i', action="append", dest='includes', default=[], type=os.path.realpath,
                        help='''base paths to search for cmd/ prep/ and parm/ directories and files''')
    parser.add_argument('-p', dest='parameters', action="append", default=[], type=os.path.realpath,
                        help='''path to additional parm/frcmod file to convert''')
    parser.add_argument('-t', action='append', dest='templates', default=[], type=os.path.realpath,
                        help='''path to additional template files to convert''')
    parser.add_argument('-s', '--skip-res', action='append', dest='skipres', default=[ "HOH", "CIO", "CIM", "CIP", "IB"],
                        help='''skip residue with this name''')
    parser.add_argument('--only', type=str, action='append', dest='only', default=[],
                        help='''only keep these templates''')
    parser.add_argument('-m', type=str, dest='fftype', default=None, help='''fftype mode for templates''')
    parser.add_argument('--skip-solvents', action='store_true', dest='skipsolvents', default=False,
                        help='''skip solvents lib''')

    if(use_extra_nbequiv_files):
        parser.add_argument('-e', action='append', dest='nbequiv', default=[], type=os.path.realpath,
                            help='''amber nb equiv files''')
    parser.add_argument('-x', action='store_true', dest='skipmissing', default=False,
                        help='''skip missing prep files (templates will be incomplete)''')
    parser.add_argument('--cleanTemplates', action='store_true', default=False,
                        help='''remove templates with missing atomtypes from templates file''')

    print("CMD: "+" ".join(sys.argv))
    args = parser.parse_args()
    if(args.rcpath is None and len(args.parameters)==0 and len(args.templates)==0):
        parser.error('must provide -r, -p or -t options')

    if(use_extra_nbequiv_files):
        nbfiles = args.nbequiv
    else:
        nbfiles = []

    ffAmber = setAmberForcefieldFiles(args.rcpath, args.includes, args.parameters, args.templates,
                                      args.fftype, nbfiles, args.skipmissing, args.skipsolvents)
    ffViparr = convertAmberForcefield(ffAmber, args.skipres, args.only, args.cleanTemplates)

    dirname=os.path.abspath(args.ffdir)
    if os.path.isdir(dirname):
        import shutil
        shutil.rmtree(dirname)

    viparr.ExportForcefield(ffViparr, dirname)
