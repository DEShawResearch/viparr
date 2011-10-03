from __future__ import print_function

from viparr import rigidify 
import msys
import pytest
import numpy as np

def convert(mol):
    return rigidify.merge_vsites_with_rigid_constraints(mol)

@pytest.fixture
def tip4p():
    return msys.Load('test/dms/tip4p.dms')

@pytest.fixture
def tip5p():
    return msys.Load('test/dms/tip5p.dms')

@pytest.fixture
def rigid():
    return msys.Load('test/dms/rigid.dms')

def test_tip4p(tip4p):
    new = convert(tip4p)
    re4 = new.table('rigid_explicit4')
    assert re4.nterms == new.nresidues
    assert re4.params.nparams == 1
    param = re4.params.param(0)
    geom = [param['%s%d' % (c,i)] for i in range(4) for c in 'xyz']
    geom = np.array(geom).reshape((-1,3))
    angle = msys.CalcAngle(geom[1], geom[0], geom[2]) * 180/np.pi
    assert np.isclose(angle, 104.52)

def test_tip5p(tip5p):
    new = convert(tip5p)
    re4 = new.table('rigid_explicit5')
    assert re4.nterms == new.nresidues
    assert re4.params.nparams == 1


def test_lc3_out_of_order():
    mol = msys.FromSmilesString('O')

    # constrain the water
    cons = mol.addTableFromSchema('constraint_hoh')
    param = cons.params.addParam()
    param['theta'] = 109
    param['r1'] = 1.0
    param['r2'] = 1.0
    cons.addTerm(mol.atoms, param)

    # add an lc3 vsite
    atoms = mol.atoms
    vsite = atoms[0].residue.addAtom()
    mol.atom(0).addBond(vsite)
    virt = mol.addTableFromSchema('virtual_lc3')
    param = virt.params.addParam()
    param['c1'] = 0.128
    param['c2'] = 0.128
    virt.addTerm([vsite]+atoms, param)

    wat2 = mol.clone()
    wat2.append(mol.clone([3,2,1,0]))

    new = convert(wat2)
    re4 = new.table('rigid_explicit4')
    assert re4.nterms == new.nresidues
    assert re4.params.nparams == 1

def test_rigid(rigid):
    new = convert(rigid)

