from __future__ import print_function

import sys
import msys
import viparr
import numpy as np

def compareRules(rule0,rule1):
    def check(name,v0,v1):
        if(v0 != v1): print("NOT EQUAL: ",name,v0,v1)

    check("vdw_func", rule0.vdw_func,rule1.vdw_func)
    check("vdw_comb_rule", rule0.vdw_comb_rule,rule1.vdw_comb_rule)
    check("plugins", sorted(rule0.plugins),sorted(rule1.plugins))
    check("nbfix_identifier", rule0.nbfix_identifier,rule1.nbfix_identifier)
    check("exclusionDepth", rule0.exclusions,rule1.exclusions)
    for i in range(2,min(rule0.exclusions,rule1.exclusions)+1):
        check("es_scale_1-%d"%(i),rule0.es_scale(i),rule1.es_scale(i))
        check("lj_scale_1-%d"%(i),rule0.lj_scale(i),rule1.lj_scale(i))

def checkTemplateExtras(tmpl0,atoms0,tmpl1,atoms1):
    for extra in ["exclusions", "nonPseudoBonds", "pseudoBonds", "angles", "dihedrals", "impropers", "cmaps"]:
       items0=getattr(tmpl0,extra)
       items1=getattr(tmpl1,extra)
       if(len(items0)!=len(items1)):
           print("Different number of template features found in ff0:%s ff1:%s for %s: %d != %d"%(tmpl0.name,tmpl1.name,extra,len(items0),len(items1)))
       if(len(items0)==0): continue
       
def compareTemplates(_templates0,_templates1):
    templates0=dict([(t.name,t) for t in _templates0])
    templates1=dict([(t.name,t) for t in _templates1])

    keys0=set(templates0.keys())
    keys1=set(templates1.keys())
    dmiss=keys1-keys0
    if(len(dmiss)): print("Missing Templates from ff0: ", dmiss)
    d=keys0-keys1
    if(len(d)): print("Missing Templates from ff1: ", d)

    typer=viparr.TemplateTyper()
    for k,v in templates0.items():
        typer.addTemplate(v)

    for name,tsys in templates1.items():
        resatoms=tsys.system.atoms
        matches,err=typer.matchFragment(tsys,resatoms)
        if(len(matches)>0):
            tmpl,resatoms=matches[0]
            matoms=tmpl.system.atoms
            for a0,a1 in zip(resatoms,matoms):
                if(a1.atomic_number<1): continue
                at0=tsys.btype(a0)
                at1=tmpl.btype(a1)
                if(not (np.isclose(a0.charge,a1.charge) or at0!=at1)):
                    # print "%d template(s) found for %s: %s"%(len(matches),tsys.name, str([t[0].name for t in matches]))
                    print("NOT EQUAL: template",tsys.name,(a0.name,a0.charge,at0),(a1.name,a1.charge,at1))          
            checkTemplateExtras(tsys,resatoms,tmpl,matoms)

        elif tsys.name not in dmiss:
            print("No template found in ff0 for %s (the template probably exists with different connectivity)"%(tsys.name))
        
def compareParams(tableName,params0,params1):
    from collections import defaultdict
    reversable=set(["stretch_harm","angle_harm","ureybradley_harm","dihedral_trig","vdw2"])
    
    store0=defaultdict(list)
    store1=defaultdict(list)
    for store,param in [(store0,params0),(store1,params1)]:
        for p in param:
            keys=list(p.keys())
            tname=""
            if("cmap_" not in tableName):
                assert("type" in keys)
                tname=p["type"]
            if(tableName in reversable):
                tnameAlt=" ".join(reversed(tname.split()))
                if tnameAlt < tname: tname=tnameAlt
            data={}
            for k in keys:
                if k in ["type","memo"]: continue
                data[k]=p[k]
            store[tname].append(data)
            
    t0=set(store0.keys())
    t1=set(store1.keys())
    d=t1-t0
    if(len(d)): print("Missing Parameters from table %s of ff0: %s"%(tableName, str(d)))
    d=t0-t1
    if(len(d)): print("Missing Parameters from table %s of ff1: %s"%(tableName, str(d)))

    for param in t0&t1:
        entry0=store0[param]
        entry1=store1[param]
        equal=len(entry0)==len(entry1)
        for e0,e1 in zip(entry0,entry1):
            assert(list(e0.keys())==list(e1.keys()))
            for k in e0:
                v0,v1=e0[k],e1[k]
                if(isinstance(v0,float)):
                    equal &= np.isclose(v0,v1)
                else:
                    equal &= v0==v1
        if not equal:
            if("cmap_" not in tableName):
                print("NOT EQUAL: ", tableName,param,store0[param]," != ",store1[param])
            else:
                print("NOT EQUAL: ", tableName,param,len(entry0),len(entry1))
                
def compareParamTables(ff0,ff1):
    t0=set(ff0.paramTables)
    t1=set(ff1.paramTables)

    d=t1-t0
    if(len(d)): print("Missing ParamTables from ff0: ", d)
    d=t0-t1
    if(len(d)): print("Missing ParamTables from ff1: ", d)

    for table in t0&t1:
        compareParams(table,ff0.params(table),ff1.params(table))


def compareCmapTables(ff0,ff1):
    tables0=ff0.cmap_tables
    tables1=ff1.cmap_tables

    if(len(tables0) != len(tables1)): print("CmapTables are of different lengths")
    for i,(t0,t1) in enumerate(zip(tables0,tables1)):
        compareParams("cmap_%d"%(i),t0.params,t1.params)


def compareForcefields(ff0,ff1):
    compareRules(ff0.rules,ff1.rules)
 
    if( isinstance(ff0.typer,viparr.TemplateTyper) and isinstance(ff1.typer,viparr.TemplateTyper)):
        compareTemplates(ff0.typer.templates,ff1.typer.templates)
        
    compareParamTables(ff0,ff1)

    compareCmapTables(ff0,ff1)

def parser():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('ff1', metavar='FFDIR')
    parser.add_argument('ff2', metavar='FFDIR')
    return parser

def main():
    args = parser().parse_args()
    ff0=viparr.ImportForcefield(args.ff1)
    ff1=viparr.ImportForcefield(args.ff2)
    compareForcefields(ff0,ff1)

