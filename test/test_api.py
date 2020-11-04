
import viparr
import msys
import pytest

@pytest.fixture
def E():
    return msys.Load('test/dms/e.dms')

@pytest.fixture
def Z():
    return msys.Load('test/dms/z.dms')

def test_template_typer():
    r=viparr.Rules()                        
    t=viparr.TemplateTyper()
    f=viparr.Forcefield(r,t)
    f.typer = t

def test_apply_ff_bond_stereo(E, Z):
    # inchis are different...
    e_inchi = msys.InChI(E, SNon=False)
    z_inchi = msys.InChI(Z, SNon=False)
    assert e_inchi != z_inchi
    
    # but a graph match exists
    e_graph = msys.Graph(E)
    z_graph = msys.Graph(Z)
    assert e_graph.match(z_graph)
    assert z_graph.match(e_graph)

    # we should not get a forcefield match between
    # ligand with different bond stereochemistry
    with pytest.raises(ValueError):
        viparr.ApplyLigandForcefields(E, [Z], selection='all')
    with pytest.raises(ValueError):
        viparr.ApplyLigandForcefields(Z, [E], selection='all')

    # but we do self-match
    viparr.ApplyLigandForcefields(E, [E], selection='all')
    viparr.ApplyLigandForcefields(Z, [Z], selection='all')


def test_apply_ff_match(E, Z):
    import numpy as np
    qE = [a.charge for a in E.atoms]
    qZ = [a.charge for a in Z.atoms]
    zE = [0 for a in E.atoms]
    zZ = [0 for a in Z.atoms]
    assert qE != qZ

    mol = E.clone()
    mol.append(Z)
    for a in mol.atoms: a.charge = 0

    matched = viparr.ApplyLigandForcefields(mol, [E,Z], selection='fragid 0')
    assert matched.natoms == mol.natoms
    assert [a.charge for a in matched.atoms] == qE + zZ

    matched = viparr.ApplyLigandForcefields(mol, [E,Z], selection='fragid 1')
    assert matched.natoms == mol.natoms
    assert [a.charge for a in matched.atoms] == zE + qZ

    matched = viparr.ApplyLigandForcefields(mol, [E,Z], selection='all')
    assert matched.natoms == mol.natoms
    assert [a.charge for a in matched.atoms] == qE + qZ

    def reorder(q0, q1, order):
        return list( np.array(q0 + q1)[order] )

    shuffle_order = list(range(0, mol.natoms, 2)) + list(range(1, mol.natoms, 2))
    shuffled = mol.clone(shuffle_order)
    # new indexing for how the original E/Z charges will be ordered after ApplyLigandForcefields
    f0, f1 = shuffled.updateFragids()
    forder = [ a.id for a in f0 + f1]
    match_order = np.array(shuffle_order)[forder]

    matched = viparr.ApplyLigandForcefields(shuffled, [E,Z], selection='fragid 0')
    assert matched.natoms == mol.natoms
    assert [a.charge for a in matched.atoms] == reorder(qE, zZ, match_order)

    matched = viparr.ApplyLigandForcefields(shuffled, [E,Z], selection='fragid 1')
    assert matched.natoms == mol.natoms
    assert [a.charge for a in matched.atoms] == reorder(zE, qZ, match_order)

    matched = viparr.ApplyLigandForcefields(shuffled, [E,Z], selection='all')
    assert matched.natoms == mol.natoms
    assert [a.charge for a in matched.atoms] == reorder(qE, qZ, match_order)

def test_apply_ff_vsite(E):
    mol = E.clone()
    mol.append(E)
    for a in mol.atoms: a.charge=0
    cl = E.select('atomicnumber 17')[0]
    v = cl.residue.addAtom()
    v.charge = 0.1
    cl.charge -= 0.1
    cl.addBond(v).order=0
    qE = [a.charge for a in E.atoms]

    matched = viparr.ApplyLigandForcefields(mol, [E])
    q = [a.charge for a in matched.atoms]
    assert q == qE + qE
    
    matched = viparr.ApplyLigandForcefields(mol, [E], selection='fragid 0')
    q = [a.charge for a in matched.atoms]
    assert q[:len(qE)] == qE

    matched = viparr.ApplyLigandForcefields(mol, [E], selection='fragid 1')
    q = [a.charge for a in matched.atoms]
    assert q[-len(qE):] == qE
    

def test_apply_ff_chiral():
    mol1 = msys.Load('test/dms/chiral1a.sdf')
    mol2 = msys.Load('test/dms/chiral1b.sdf')
    s1 = msys.InChI(mol1, SNon=False).string
    s2 = msys.InChI(mol2, SNon=False).string
    assert s1 != s2
    s1 = msys.InChI(mol1, SNon=True).string
    s2 = msys.InChI(mol2, SNon=True).string
    assert s1 == s2

    # this should work because we ignore tet stereo by default
    viparr.ApplyLigandForcefields(mol1, [mol2])

    # but this makes it fail
    with pytest.raises(ValueError):
        viparr.ApplyLigandForcefields(mol1, [mol2], match_tet_stereo=True)


def testFixProchiralProteinAtomNames():
    # These are 4 NMR structures (so they have resolved protons). The PDB follows
    # the IUPAC conventions, so that means we shouldn't get any flips/
    pdb_codes = ["1D3Z", "1G03", "1XPW", "2KOD"]
    dmses = []
    for code in pdb_codes:
        dmses.append(msys.Load(f"test/dms/prochiral/{code}.0.pdb"))

    for dms in dmses:
        flippedAtomIds = viparr.FixProchiralProteinAtomNames(dms)
        assert len(flippedAtomIds) == 0

def testViparrRenameProchiral_1():
    mol = msys.Load("test/dms/ww.dms", structure_only=True)
    ffs = [viparr.ImportForcefield(viparr.find_forcefield(f))
            for f in ("aa.amber.ff99", "water.tip3p")]
    mol1 = mol.clone()
    mol2 = mol.clone()
    viparr.ExecuteViparr(mol1, ffs, rename_atoms=True, rename_residues=True, rename_prochiral_atoms=True, verbose=True)
    viparr.ExecuteViparr(mol2, ffs, rename_atoms=True, rename_residues=True, rename_prochiral_atoms=False, verbose=True)
    assert mol1.natoms==mol2.natoms
    renamed = [(a,b) for a,b in zip(mol1.atoms, mol2.atoms) if a.name != b.name]
    assert len(renamed) == 14

def testViparrRenameProchiral_2():
    mol = msys.Load("test/dms/prochiral/VVV.dms", structure_only=True)
    ids = viparr.FixProchiralProteinAtomNames(mol)
    names = [f"{mol.atom(i).residue.name}:{mol.atom(i).name}" for i in ids]
    assert names == ["VAL:CB", "VAL:CB", "VAL:CB"]
