from __future__ import print_function

import msys
import numpy as np
import re

def construct_ahnr_geometry(param, N):
    n = N*(N+1)//2
    r = [param['r%d' % i] for i in range(1, n+1)]
    if len(r) > 3:
        return new_place_ahnr(*r)
    else:
        return place_ahnr(*r)

def construct_h2o_geometry(param):
    theta = param['theta'] * np.pi / 180
    r1 = param['r1']
    r2 = param['r2']
    x = np.cos(theta) * r2
    y = np.sin(theta) * r2
    geometry = list()
    geometry.append([0,0,0])
    geometry.append([r1,0,0])
    geometry.append([x,y,0])    

    return np.array(geometry)

def construct_rigid_geometry(param):
    pos = [param[k] for k in sorted(param.keys())]
    pos = np.array(pos).reshape(3,-1).T
    return pos

def construct_lcn(param, n, pos):
    c = [param['c%d' % i] for i in range(1, n)]
    return place_lcn(c, pos)


def normalized(x):
    """Return normalized copy of vector x
    """
    return np.array(x) / np.linalg.norm(x)

def place_ahnr(*r):
    """ construct geometry from ahNR params using linear least squares"""
    from scipy.optimize import leastsq
    from scipy.spatial.distance import pdist

    def pos_from(x):
        """ map free parameters to spatial coordinates """
        assert len(x) == 1 or len(x) % 3 == 0
        n = len(x) // 3 + 2
        pos = np.zeros(n*3) # p0 at the origin
        pos[3] = x[0]       # p1 aligned with x axis
        pos[6:8] = x[1:3]   # p2 in xy plane
        pos[9:] = x[3:]     # all other points
        return pos.reshape((n,3))

    def f(x, *r):
        """ penalty function """
        return pdist(pos_from(x)) - r

    x0 = np.zeros(len(r))
    x = leastsq(f, x0, r)[0]
    return pos_from(x)

def new_place_ahnr(*r):
    """ construct geometry from ahNR params using multivariate optimization"""
    from scipy.optimize import minimize
    from scipy.spatial.distance import pdist

    def f(x, *r):
        """ penalty function """
        y = pdist(x.reshape(-1, 3)) - r
        return np.dot(y, y)

    R = len(r)
    N = int((1+8*R)**0.5 + 1)//2
    x0 = np.zeros(3*N)
    opt = minimize(f, x0, r, method='Powell')
    assert opt.success, opt
    return opt.x.reshape((-1,3))

def place_lcn(c, pos):
    ''' lifted from zero/mechanics/virtuals/virtual_lc.hxx '''
    N = len(c) + 1
    assert len(pos)==N, "Expected %d positions; got %d" % (len(pos), N)
    vpos = np.zeros(3)
    csum = 1.0
    for i in range(1,N):
        vpos += c[i-1] * pos[i]
        csum -= c[i-1]
    return csum * pos[0] + vpos

def place_out3(c, pos):
    N = 3
    assert len(pos)==N, "Expected %d positions; got %d" % (len(pos), N)
    assert len(c)==N, "Expected %d coefficients; got %d" % (len(c), N)
    a = pos[1] - pos[0]
    b = pos[2] - pos[0]
    c1, c2, c3 = c
    return pos[0] + c1*a + c2*b + c3*np.cross(a,b)

def place_out3n(c, pos):
    N = 3
    assert len(pos)==N, "Expected %d positions; got %d" % (len(pos), N)
    assert len(c)==N, "Expected %d coefficients; got %d" % (len(c), N)
    a = pos[1] - pos[0]
    b = pos[2] - pos[0]
    q = np.cross(a, b)
    c1, c2, c3 = c
    return pos[0] + c1*normalized(a) + c2*normalized(b) + c3*normalized(q)

def add_rigid_explicit(N, mol):
    assert N >= 2, "N must be at least 2"
    table = mol.addTable('rigid_explicit%d' % N, N)
    table.category = 'constraint'
    params = table.params
    for i in range(N):
        params.addProp('x%d' % i, float)
        params.addProp('y%d' % i, float)
        params.addProp('z%d' % i, float)
    return table

def merge_vsites_with_rigid_constraints(mol, include_rigid=False):
    mol = mol.clone()
    vtables = [table for table in mol.tables if table.category == 'virtual']
    reich_pattern = re.compile(r'constraint_ah([1-9])R')
    lcn_pattern = re.compile(r'virtual_lc([1-9])')
    for table in mol.tables:
        reich = reich_pattern.match(table.name)
        if table.name == 'constraint_hoh':
            print("Construct geometry from", table.name)
            gmap = { p.id : construct_h2o_geometry(p) for p in table.params.params }
        elif reich:
            print("Construct geometry from", table.name)
            N = int(reich.group(1))
            gmap = { p.id : construct_ahnr_geometry(p, N) for p in table.params.params }
        elif table.name.startswith('rigid_explicit'):
            if include_rigid:
                gmap = { p.id : construct_rigid_geometry(p) for p in table.params.params }
            else:
                print("Skipping", table.name)
                continue
        elif table.category == 'constraint':
            print("Warning, skipping non-rigid constraint table '%s'" % table.name)
            continue
        else:
            continue
        print("Processing %d terms from %s" % (table.nterms, table.name))
        for term in table.terms:
            term_atoms = term.atoms
            geometry = gmap[term.param.id]
            vsites = dict()
            for vtable in vtables:
                lcn = lcn_pattern.match(vtable.name)
                for vterm in vtable.findWithAny(term_atoms):
                    vatoms = vterm.atoms
                    # vsite must be completely defined by atoms in the constraint group
                    if not set(vatoms[1:]).issubset(term_atoms):
                        continue
                    pos = [geometry[term_atoms.index(a)] for a in vatoms[1:]]
                    param = vterm.param
                    if lcn:
                        N = int(lcn.group(1))
                        vpos = construct_lcn(param, N, pos)
                    elif vtable.name == 'virtual_out3':
                        c = [param['c%d' % i] for i in range(1, len(term_atoms)+1)]
                        vpos = place_out3(c, pos)
                    elif vtable.name == 'virtual_out3n':
                        c = [param['c%d' % i] for i in range(1, len(term_atoms)+1)]
                        vpos = place_out3n(c, pos)
                    else:
                        raise RuntimeError("Unsupported vsite table '%s'" % vtable.name)
                    vsites[vatoms[0].id] = vpos
                    vterm.remove()

            N = len(geometry) + len(vsites)
            etable = add_rigid_explicit(N, mol)
            param = etable.params.addParam()
            atoms = term_atoms
            for i, (x,y,z) in enumerate(geometry):
                param['x%d' % i] = x
                param['y%d' % i] = y
                param['z%d' % i] = z
            for i, id in enumerate(sorted(vsites)):
                x,y,z = vsites[id]
                param['x%d' % len(atoms)] = x
                param['y%d' % len(atoms)] = y
                param['z%d' % len(atoms)] = z
                atoms.append(mol.atom(id))
            etable.addTerm(atoms, param)
        table.remove()

    for vtable in vtables:
        if vtable.nterms == 0:
            vtable.remove()

    mol.coalesceTables()
    new = mol.clone()

    for table in new.tables:
        if table.name.startswith('rigid_explicit'):
            print("%18s: %5d terms, %d distinct params" % (table.name, table.nterms, table.params.nparams))
    return new


