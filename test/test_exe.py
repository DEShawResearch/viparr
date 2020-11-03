import subprocess
import os
import shlex
import pytest
import msys
import numpy

here=os.path.dirname(__file__)
dmsdir=f'{here}/dms'
ffdir=f'{here}/ff3'
amberdir=f'{here}/conversion/amber'
chmdir=f'{here}/conversion/charmm'

def call(cmd):
    if isinstance(cmd, str):
        cmd = shlex.split(cmd)
    subprocess.run(cmd, check=True)

def fail(cmd):
    with pytest.raises(subprocess.CalledProcessError):
        call(cmd)


def test_iviparr(tmpdir):
    subprocess.check_call(['iviparr', '--ffname', 'aa.amber.ff99SB', '%s/e.dms' % dmsdir, str(tmpdir.join('out.ff'))])

def test_build_constraints(tmpdir):
    subprocess.run(["viparr-build-constraints", f"{dmsdir}/ww.dms", "-o", f"{tmpdir}/out.dms"], check=True)

def test_compare_forcefields(tmpdir):
    subprocess.run(["viparr-compare-forcefields", f"{ffdir}/amber03", f"{ffdir}/amber99"], check=True)

def test_convert_amber(tmpdir):
    subprocess.run(f"viparr-convert-amber -p {amberdir}/parm10.dat -p {amberdir}/frcmod.ff14SB -t {amberdir}/amino12.in -t {amberdir}/aminoct12.in -t {amberdir}/aminont12.in -m amino_acids {tmpdir}/ff".split(), check=True)
    subprocess.run(["diff", f"{tmpdir}/ff", f"{amberdir}/convert_reference"], check=True)

def test_convert_charmm(tmpdir):
    subprocess.run(f"viparr-convert-charmm -p {chmdir}/par_all36m_prot.prm -n {chmdir}/amino_acids -t {chmdir}/top_all36_prot.rtf {tmpdir}/ff".split(), check=True)
    subprocess.run(["diff", f"{tmpdir}/ff", f"{chmdir}/convert_reference"], check=True)

def test_solvate(tmpdir):
    # mismatched forcefield
    fail(f"viparr_solvate {dmsdir}/ww.dms {tmpdir}/solv.dms --ffname water.tip3p")
    # correct forcefield
    call(f"viparr_solvate {dmsdir}/ww.dms {tmpdir}/solv.dms --ffname water.spc_opls")
    # make sure coordinates weren't changed by solvate
    oldpos = msys.Load(f"{dmsdir}/ww.dms").positions
    newpos = msys.Load(f"{tmpdir}/solv.dms").positions
    assert numpy.allclose(oldpos, newpos[:len(oldpos)])

    # no forcefield ok if none in input system
    call(f"viparr_solvate {dmsdir}/chiral1a.sdf {tmpdir}/solv.dms")

    # FIXME - this should fail because ww.dms has a forcefield.  we don't want solvate to silently lose information
    # no forcefield
    #fail(f"viparr_solvate {dmsdir}/ww.dms {tmpdir}/solv.dms")
