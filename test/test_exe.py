import subprocess
import os
here=os.path.dirname(__file__)
dmsdir=f'{here}/dms'
ffdir=f'{here}/ff3'
amberdir=f'{here}/conversion/amber'
chmdir=f'{here}/conversion/charmm'

def test_iviparr(tmpdir):
    import garden
    garden.load('viparr-ff/2.1.0c7/data')
    subprocess.check_call(['iviparr', '--ffname', 'aa.amber.ff99SB', '%s/e.dms' % dmsdir,
            str(tmpdir.join('out.ff'))])

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
