import os, msys, viparr
import numpy

mydir = os.path.dirname(__file__)
WAT=os.path.join(mydir, '..', '..', '..', '..', 'share', 'solvate', 'h2o.dms')
TIP3PWAT=os.path.join(mydir, '..', '..', '..', '..', 'share', 'solvate', 'h2o.tip3p.dms')
WATRAD = 2.4
WATSEL = 'oxygen'
WATCON = 1.0

def remove_periodic_contacts(_mol, npro, dist):
    mol = _mol.clone()
    natoms = mol.natoms
    wat = _mol.clone('index >= %d' % npro)
    mol.append(wat)
    sel='index < %(natoms)d and index >= %(npro)d and within %(dist)f of (index >= %(natoms)d)' % locals()
    pos = mol.positions
    a, b, c = mol.cell
    bad = []
    for i in (-1,0,1):
        for j in (-1,0,1):
            for k in (-1,0,1):
                if i==0 and i==j and i==k: continue
                delta = i*a + j*b + k*c
                pos[natoms:] += delta
                mol.setPositions(pos)
                ids = mol.selectIds(sel)
                bad.extend(ids)
                pos[natoms:] -= delta

    if bad:
        _mol=_mol.clone('not same fragid as index ' + ' '.join(map(str, bad)))
    return _mol

def Solvate(mol, watbox=None, dims=None, center=None,
            chain='WT', verbose=True, ffname='', ffdir='', ff=[],
            without_constraints=False, water_free_zone_width = WATCON,
            without_fix_masses=True, solvent_radius=WATRAD,
            solvent_selection=WATSEL):
    ''' 
    Build and return a new system consisting of mol plus solvent.

    If no box dimensions are provided, the box size along each axis will
    be equal to the greater of that of the input structure and the input
    water box.

    Return the solvated system; no modifications are made to the input system.
    '''

    mol=mol.clone()
    npro = mol.natoms

    # put all the water in one ct
    ct=mol.addCt()
    ct.name='solvate'
    if watbox is None:
        wat=msys.Load(TIP3PWAT, structure_only=True)
    elif ffdir != '' or ffname != '':
        wat=msys.Load(watbox, structure_only=True)
    else:
        wat=msys.Load(watbox, structure_only=False)

    if verbose:
        print("Loaded water box from '%s'" % wat.name)
    watsize=[wat.cell[i][i] for i in range(3)]
    molsize=[mol.cell[i][i] for i in range(3)]

    if ffdir != '' or ffname != '' or ff:
        if ff:
            if isinstance(ff,viparr.Forcefield):
                ff = [ff]
            else:
                ff = list(ff)
        else:
            # retain backward compatible behavior of ffdir, accepting relative path
            ffpath = ffdir or viparr.find_forcefield(ffname)
            ff = [viparr.ImportForcefield(ffpath)]
            if verbose:
                print('Applying forcefield ' + ffpath + ' to water box')
        viparr.ExecuteViparr(wat, ff,
                with_constraints=(not without_constraints), verbose=verbose)

    # find a chain name for the waters that doesn't overlap with the 
    # input structure.
    chains = set(c.name for c in mol.chains)
    watchain = 'X'
    while watchain in chains:
        watchain += 'X'
    for c in wat.chains:
        c.name = watchain

    if verbose:
        print("Dims specified as '%s'" % dims)
    if dims is None:
        dims = [max(x) for x in zip(watsize, molsize)]
    elif len(dims)==1:
        dims=[float(dims[0])]*3
    elif len(dims)!=3:
        raise ValueError("Dims must be given as list of one or three values")
    else:
        dims=[float(x) for x in dims]
    
    if center is None:
        center=[0,0,0]
    elif len(center)!=3:
        raise ValueError("Center must be given as list of three values")
    else:
        center=[float(x) for x in center]

    # update the cell
    mol.cell[0][:]=[dims[0],0,0]
    mol.cell[1][:]=[0,dims[1],0]
    mol.cell[2][:]=[0,0,dims[2]]

    # extract where to put the water
    xmin = center[0]-0.5*dims[0]
    ymin = center[1]-0.5*dims[1]
    zmin = center[2]-0.5*dims[2]
    xmax = center[0]+0.5*dims[0]
    ymax = center[1]+0.5*dims[1]
    zmax = center[2]+0.5*dims[2]
    nx = int(dims[0]/watsize[0]) + 1
    ny = int(dims[1]/watsize[1]) + 1
    nz = int(dims[2]/watsize[2]) + 1

    xshift = -0.5 * (nx-1)*watsize[0]
    yshift = -0.5 * (ny-1)*watsize[1]
    zshift = -0.5 * (nz-1)*watsize[2]

    # replicate the template water box
    if verbose: print("replicating %d x %d x %d" % (nx,ny,nz))
    for i in range(nx):
        xdelta = xshift + i*watsize[0]
        for j in range(ny):
            ydelta = yshift + j*watsize[1]
            for k in range(nz):
                zdelta = zshift + k*watsize[2]
                newatoms = ct.append(wat)
                for a in newatoms:
                    a.x += xdelta
                    a.y += ydelta
                    a.z += zdelta
    mol.updateFragids()

    if verbose: print("removing overlaps")

    toonear='pbwithin %s of index < %s' % (solvent_radius, npro)

    # remove overlap with solute
    mol=mol.clone(
        'not same fragid as (index >= %d and (%s) and (%s))' % (
            npro, solvent_selection, toonear))
    if verbose: print("After removing overlap, %d solvent atoms" % (
            mol.natoms - npro))

    # remove molecules whose center is outside the desired box
    hits = mol.select('index >= %d and (x<%s or y<%s or z<%s or x>%s or y>%s or z>%s)' % (npro, xmin,ymin,zmin, xmax,ymax,zmax))
    frags = set(a.fragid for a in hits)
    if frags:
        fragmap = dict()
        sel = '(%s) and fragid %s' % (solvent_selection, ' '.join(map(str,frags)))
        for a in mol.select(sel): fragmap.setdefault(a.fragid, []).append(a.id)
        outside = []
        pos = mol.getPositions()
        pmin = numpy.array((xmin,ymin,zmin))
        pmax = numpy.array((xmax,ymax,zmax))
        for fragid, ids in fragmap.items():
            c = pos[ids].mean(0)
            if (c<pmin).any() or (c>pmax).any(): outside.append(fragid)
        if outside:
            mol = mol.clone('not fragid ' + ' '.join(map(str,outside)))
            if verbose: print("After removing outside solvent molecules, %d solvent atoms" % (mol.natoms - npro))

    # remove overlap with periodic images
    mol = remove_periodic_contacts(mol, npro, water_free_zone_width)
    if verbose: print("after removing periodic clashes, %d solvent atoms" % (mol.natoms-npro))

    # assign the water chain name and water resids
    watres = 1
    for c in mol.chains:
        if c.name == watchain:
            c.name = chain
            for r in c.residues:
                r.resid = watres
                watres += 1

    # fix mass
    if not without_fix_masses:
        viparr.FixMasses(mol, mol.select('water'), verbose=verbose)

    if verbose: print("updating global cell to (%g %g %g)" % tuple(dims))

    return mol

__doc = \
'''viparr_solvate watbox_out.dms [ options ]  -- create water box
viparr_solvate solute_in.dms solvate_out.dms [ options ] -- add water to solute
viparr_solvate solute_in.dms watbox_in.dms solvate_out.dms [ options ] -- add
    specified water box to solute

Generate water molecules around a structure. With one argument, a system
containing only water is created. With two arguments, a generic water box
is tiled around the input structure. With three arguments, the provided
water box is tiled around the input structure.

If a forcefield is specified, the generic or provided water box is
parameterized with the forcefield. Otherwise, the generic water box is
unparameterized and the provided water box keeps its original parameters.

If thickness T is specified, then a cubic water box is generated with
dimensions equal to W + 2*T, where W is the greatest extent of the
solute along any coordinate axis.  Only one of dims and thickness may
be specified.

If minimization is specified, the solute will be rotated such that the
cubic box circumscribing it has minimum volume. This option can be
used with the thickness option. It is an error to specify this option
with dims.

By default, waters within 1A of a periodic neighbor will be removed.
This should not produce significant vacuum regions in your system;
you should, however, always equilibrate your system after solvating it.

'''

import sys, os
import argparse

class RawDescriptionArgumentDefaultsFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def main():
    parser = argparse.ArgumentParser(description=__doc, formatter_class=RawDescriptionArgumentDefaultsFormatter)

    parser.add_argument('dms', nargs='+',
            help="Between one and three dms files. See above for how these are interpreted.")
    parser.add_argument('-d', '--dims', default=None,
            help='water box dimensions: 1 or 3 comma-separated values')
    parser.add_argument('-t', '--thickness', default=None, type=float,
            help='Minimum distance between solute and edge of water box')
    parser.add_argument('-c', '--center', default=None,
            help='center of box as 3 comma-separated values; default 0,0,0')
    parser.add_argument('-n', '--chain', default='WT',
            help='Chain name of constructed water box')
    parser.add_argument('-v', '--verbose', default=False, action='store_true')
    parser.add_argument('--ffname', default='',
            help='forcefield name, if forcefield is to be applied to water box')
    parser.add_argument('--ffdir', default='', help="""explicit forcefield
    directory, if forcefield is to be applied to water box""")
    parser.add_argument('--without-constraints', default=False,
           action='store_true', help='apply forcefield without constraints')
    parser.add_argument('-w', '--water-free-zone-width', default=1.0, type=float,
            help='Minimum distance between any water molecule and the edge of the box')
    parser.add_argument('--without-fix-masses', action='store_true', default=False,
            help='Disable assignment of identical mass to all waters')
    parser.add_argument('--solvent-radius', type=float, default=WATRAD,
            help='Radius of solvent molecules used in solvent-solute distance check')
    parser.add_argument('--solvent-selection', default=WATSEL,
            help='Selection of solvent atoms used in solvent-solute distance check')

    opts = parser.parse_args()
    args = opts.dms[:]
    del opts.dms
    watbox = None

    if opts.ffname != '' and opts.ffdir != '':
        parser.error("must specify at most one forcefield")
    if opts.dims is not None and opts.thickness is not None:
        parser.error("cannot specify both dims and thickness")
    if len(args)==1:
        mol=msys.CreateSystem()
        output=args[0]
    elif len(args)==2:
        input, output = args
        if opts.verbose: print("Loading input file <%s>" % input)
        if opts.ffname != '' or opts.ffdir != '':
            mol=msys.LoadDMS(input, structure_only=False)
        else:
            mol=msys.LoadDMS(input, structure_only=True)
    elif len(args)==3:
        input, watbox, output = args
        if opts.verbose: print("Loading input file <%s>" % input)
        mol=msys.LoadDMS(input, structure_only=False)
    else:
        parser.error("incorrect number of arguments")

    if opts.dims is not None:
        opts.dims = [float(x) for x in opts.dims.split(',')]
    if opts.thickness is not None:
        pos = mol.getPositions()
        extent = max(pos.max(0)-pos.min(0))
        extent += 2*opts.thickness
        extent = float(extent)
        opts.dims = [extent,extent,extent]
        print("Final box dimensions:", opts.dims)
    del opts.__dict__['thickness']
    if opts.center is not None:
        opts.center = [float(x) for x in opts.center.split(',')]

    mol = Solvate(mol, watbox, **opts.__dict__)

    # create nonredundant parameter tables
    mol.coalesceTables()
    mol=mol.clone()

    if opts.verbose: print("Writing DMS file <%s>" % output)
    msys.SaveDMS(mol,output)

