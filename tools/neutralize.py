"""
Adds or removes ions from a chemical system in order to neutralize the charge
and attain a specified ionic concentration. A forcefield may be provided to
parametrize any added ions. Randomly selected water molecules may be removed
from the system to make room for the added ions. In the following, counterions
are those with opposite charge from the solute (required to neutralize the
charge of the system), and counter-counterions are those of the same charge as
the solute (possibly required to achieve a desired ionic concentration).
"""

import os, msys, viparr, math
from collections import namedtuple
import random

IONPROPS = {
        # FIXME: we should get the list of recognized ion species
        # from the input forcefield.  until then, add a few more
        # names.
        'NA' :  (11,1),
        'Na+':  (11,1),
        'K'  :  (19,1),
        'K+' :  (19,1),
        'CL' :  (17,-1),
        'Cl-':  (17,-1),
        'MG' :  (12,2),
        'Mg2+': (12,2),
        'CA' :  (20,2),
        'Ca2+': (20,2),
        }

def parse_ff(ffname, ffdir, ff):
    if ffdir != '' or ffname != '' or ff:
        if ff:
            if isinstance(ff,viparr.Forcefield):
                ff = [ff]
            else:
                ff = list(ff)
        else:
            ffpath = ffdir or viparr.find_forcefield(ffname)
            ff = viparr.ImportForcefield(ffpath)
    return ff

def parse_ion(name, ff, verbose):
    ionprop = IONPROPS.get(name)
    if ionprop is None: ionprop = IONPROPS.get(name.upper())
    if ionprop is None:
        raise ValueError("Unrecognized ion type '%s'" % name)
    anum, charge = ionprop
    ionsys = msys.CreateSystem()
    chain = ionsys.addChain()
    res = chain.addResidue()
    res.name = name
    atom = res.addAtom()
    atom.name = name
    atom.atomic_number=anum
    atom.charge=charge
    atom.mass = msys.MassForElement(anum)
    if ff:
        viparr.ExecuteViparr(ionsys, [ff], verbose=verbose)

    return ionsys

def compute_center(residue):
    tm=0.0
    tx=0.0
    ty=0.0
    tz=0.0
    for a in residue.atoms:
        m = a.mass
        tm += m
        tx += m*a.x
        ty += m*a.y
        tz += m*a.z
    if tm:
        tx /= tm
        ty /= tm
        tz /= tm
    return (tx,ty,tz)

def dist2(pi, pj):
    d=0.0
    for i,j in zip(pi, pj):
        d += (i-j)**2
    return math.sqrt(d)


def Neutralize(mol, cation='NA', anion='CL', 
        chain='ION', chain2='ION2',
        solute_pad=5.0,
        ion_pad=3.0,
        concentration=0.0,
        keep='none',
        verbose=True,
        ffname='',
        ffdir='',
        ff=[],
        random_seed=0):

    """
    neutralize() -> replace water molecules with ions

    Will build and return a new system, keeping the input system the same.
  
    Optional arguments (default):
      cation  ('NA') -- species of cation.  Must be NA or K
      anion   ('CL') -- species of anion.  Currently only CL is supported
      chain ('ION'), chain2 ('ION2') -- counterions (those with the opposite
                          charge of the solute) are placed in segment ION; the
                          other species is put in segment ION2.
      solute_pad (5) -- minimum distance between placed ions and non-water.
      ion_pad    (3) -- minimum distance between placed ions.
      concentration (0) -- molar concentration of added ions with same charge
                          as the solute.  The number of such ions will be
                          given by 

                            nother = int((conc / 55.345) * (ntotalwat - nions))

                          where ntotalwat is the total number of waters in
                          the original system and nions is the number of 
                          counterions needed to achieve charge neutrality.
      keep  ('none') -- atomsel of ions/waters that cannot be deleted or replaced
      random_seed (0) -- seed random number generator to determine which water
                       molecules to replace

    """

    # Clone system
    mol = mol.clone()

    # create a single ct for all the ions
    ct = mol.addCt()
    ct.name = 'ion'

    # first, we add sufficient counterions to neutralize.  Then we add 
    # counterions and counter-counterions until counter-counterions are
    # up to the desired concentration.  
    ff = parse_ff(ffname, ffdir, ff)
    cationsys = parse_ion(cation, ff, verbose)
    anionsys = parse_ion(anion, ff, verbose)
    cationatom = cationsys.atoms[0]
    anionatom = anionsys.atoms[0]

    solute=mol.select("""not ((atomicnumber %d and not bonded) or (atomicnumber %d and not bonded)) or (%s)""" %
                      (cationatom.atomic_number, anionatom.atomic_number, keep))
    cg = sum(a.charge for a in solute)
    if verbose:
        print("neutralize: solute charge=%s" % cg)

    if cg >= 0:
        ionsys = anionsys
        othersys = cationsys
        ionatom = anionatom
        otheratom = cationatom
    else:
        ionsys = cationsys
        othersys = anionsys
        ionatom = cationatom
        otheratom = anionatom

    nions = int(math.fabs(cg/ionatom.charge)+0.5)

    # find the water residues
    water = mol.select('water and (not hydrogen) and (not within %f of (not water)) and (not (%s))'
            % (solute_pad, keep))
    residues = sorted(set(a.residue for a in water), key=lambda x: x.id)
    nwat = len(residues)
    if verbose:
        print("waters available to be replaced by ions:", nwat)

    # compute number of ions already present
    nions_prev = len(mol.select('atomicnumber %d and not bonded and (not (%s))'
                                % (ionatom.atomic_number, keep)))
    nother_prev = len(mol.select('atomicnumber %d and not bonded and (not (%s))'
                                 % (otheratom.atomic_number,keep)))

    if verbose:
        print("Starting with %d %s ions" % (nions_prev, ionatom.name))
        print("Starting with %d %s ions" % (nother_prev, otheratom.name))

    # convert molar concentration to number based on available waters.  The
    # molar concentration of water is about 55.345 mol/L.  Use all available
    # waters to calculate the number of ions to add.
    ntotalwat = len(set(a.residue for a in mol.select('water')))
    if verbose:
        print("Starting with %d water molecules" % ntotalwat)
    cgratio = math.fabs(otheratom.charge/ionatom.charge)
    nother = int((concentration/55.345) * (ntotalwat-nions+nions_prev))
    nions += int(nother*cgratio)

    # subtract off the ions already present in solution
    nions -= nions_prev
    nother -= nother_prev

    qtot = sum(a.charge for a in mol.atoms)
    qnew = qtot + nother*otheratom.charge + nions*ionatom.charge
    if abs(qnew) > 0.001:
        # this can happen if cgratio > 1.
        if abs(ionatom.charge) < abs(otheratom.charge):
            nions -= round(qnew / ionatom.charge)
        else:
            nother -= round(qnew/otheratom.charge)

    if nions >= 0 and verbose:
        print("New system should contain %d %s ions" % (nions, ionatom.name))
    if nother >= 0 and verbose:
        print("New system should contain %d %s ions" % (nother, otheratom.name))

    if nions < 0:
        # delete ions
        sel=mol.select('atomicnumber %d and not bonded and not (%s)' %
                (ionatom.atomic_number, keep))
        if len(sel) < -nions:
            raise RuntimeError("""Cannot decrease concentration - 
            not enough ions to delete""")
        for r in sel[:-nions]: r.remove()
        nions = 0
    if nother < 0:
        # delete other 
        sel=mol.select('atomicnumber %d and not bonded and not (%s)' %
                (otheratom.atomic_number, keep))
        if len(sel) < -nother:
            raise RuntimeError("""Cannot decrease concentration - not enough 
            other ions to delete""")
        for r in sel[:-nother]: r.remove()
        nother = 0


    if nwat < nions + nother:
        raise RuntimeError("""Only %d waters found; not enough to 
        neutralize""" % nwat)

    # Shuffle the residues
    random.seed(random_seed)
    random.shuffle(residues)

    # Remove overly close residues among the first nions + nother waters
    # Cache the centers to save time
    centers={}
    ionpad2 = ion_pad * ion_pad
    for i in range(nions + nother):
        ri = residues[i]
        try: pi = centers[ri.id]
        except KeyError:
            pi = centers[ri.id] = compute_center(ri)
        j = i+1
        while j < nions+nother:
            if len(residues) < nions+nother:
                raise RuntimeError("Not enough waters or too large ion_pad.")
            rj = residues[j]
            try: pj = centers[rj.id]
            except KeyError:
                pj = centers[rj] = compute_center(rj)
            d2 = dist2(pi,pj)
            if d2 < ionpad2:
                del residues[j]
            else:
                j += 1

    keep_atoms = set([x.id for x in mol.select('(%s)' % (keep,))])
    keep_residues = set([x.residue.id for x in mol.select('(%s)' % (keep,))])
    residues_removed = set()

    def multiple_waters(residue):
        real_atoms = [atom for atom in residue.atoms if atom.atomic_number > 0]
        return len(real_atoms) > 3

    def handle_multiple_waters(residue):
        raise RuntimeError("Some residues (%d) contain multiple waters. Use dms-fix-water-residues." % residue.id)

    if nions > 0:
        for i in range(nions):
            res=residues[i]
            if multiple_waters(res):
                handle_multiple_waters(res)
            newion = ct.append(ionsys)[0]
            newion.pos = compute_center(res)
            newion.residue.chain.name = chain
            newion.residue.resid=i+1
            residues_removed.add(res.id)
            res.remove()

    if nother > 0:
        for i in range(nions, nions+nother):
            res=residues[i]
            if multiple_waters(res):
                handle_multiple_waters(res)
            newion = ct.append(othersys)[0]
            newion.pos = compute_center(res)
            newion.residue.chain.name = chain2
            newion.residue.resid=i+1-nions
            residues_removed.add(res.id)
            res.remove()

    assert not (residues_removed & keep_residues)

    if ff and ff.rules.nbfix_identifier:
        if verbose:
            print("nbfix terms in ion forcefield; re-running viparr to generate overrides")
        ffs = [ff]
        fftable = mol.auxtable('forcefield')
        for p in fftable.params:
            ffs.append(viparr.ImportForcefield(p['path']))
        for table in mol.tables: table.remove()
        mol = mol.clone('atomicnumber > 0')
        viparr.ExecuteViparr(mol, ffs, verbose=verbose, verbose_matching=False)
    else:
        mol.coalesceTables()

    qtot = sum(a.charge for a in mol.atoms)
    assert abs(qtot) < 0.01, f"Neutralization failed, total charge after neutralization is {qtot:.3f}!"
    return mol.clone()

__doc = \
'''
viparr_neutralize input.dms output.dms [ options ]

Add/remove counterions to a structure in order to neutralize the charge and
attain a specified ionic concentration. A forcefield may be provided to
parameterize the added ions.
'''

IONS = IONPROPS
cations = ', '.join(x for x in IONS if IONS[x][1] > 0)
anions  = ', '.join(x for x in IONS if IONS[x][1] < 0)

def main():
    import argparse
    parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('ifile', help='input DMS file')
    parser.add_argument('ofile', help='output DMS file')
    parser.add_argument('-p', '--cation', default='NA', 
            help='Species of cation: one of %s' % cations)
    parser.add_argument('-n', '--anion', default='CL',
            help='Species of anion: one of %s' % anions)
    parser.add_argument('-c', '--chain', default='ION',
            help='Chain name for counterions')
    parser.add_argument('-C', '--chain2', default='ION2',
            help='Chain name for counter-counterions')
    parser.add_argument('-s', '--solute-pad', default=5.0, type=float,
            help='minimum distance between placed ions and solute')
    parser.add_argument('-i', '--ion-pad', default=3.0, type=float,
            help='minimum distance between placed ions')
    parser.add_argument('-m', '--concentration', default=0.0, type=float,
            help='molar concentration of counter-counterions')
    parser.add_argument('-k', '--keep', default='none',
            help='Atomsel of ions/waters that should not be deleted or replaced')
    parser.add_argument('-v', '--verbose', default=False, action='store_true')
    parser.add_argument('--no-ff', default=False, action='store_true',
            help='Delete all forcefield information')
    parser.add_argument('--ffname', default='',
            help='forcefield name, if forcefield is to be applied to new ions')
    parser.add_argument('--ffdir', default='',
            help="""explicit forcefield directory, if forcefield is to be 
            applied to new ions""")
    parser.add_argument('--random-seed', default=0,
            help="""seed for random number generator to determine which
            water molecules to replace""")

    args = parser.parse_args()

    if args.no_ff:
        if args.ffname or args.ffdir:
            parser.error("Cannot specify both --no-ff and a forcefield")
    else:
        if bool(args.ffname) == bool(args.ffdir):
            parser.error("Must specify exactly one forcefield, or --no-ff")

    if args.verbose:
        print("Loading input file <%s>" % args.ifile)
    mol = msys.Load(args.ifile, structure_only=args.no_ff)
    cfg = args.__dict__
    del cfg['no_ff']

    if not mol.table_names:
        parser.error("\nERROR: This system has not yet been parameterized with a \
force-field. You must run viparr_neutralize AFTER \
applying your force-field, because force-field application may \
change the net charge.")

    mol = Neutralize(mol, **cfg)
        
    # create nonredundant parameter tables
    mol.coalesceTables()
    mol=mol.clone()

    if args.verbose:
        print("Writing DMS file <%s>" % args.ofile)
    msys.SaveDMS(mol, args.ofile)

