from __future__ import print_function

import subprocess
import os
dmsdir='%s/dms' % os.path.dirname(__file__)
ffdir='%s/ff3' % os.path.dirname(__file__)

def test_iviparr(tmpdir):
    import garden
    garden.load('viparr-ff/2.1.0c7/data')
    subprocess.check_call(['iviparr', '--ffname', 'aa.amber.ff99SB', '%s/e.dms' % dmsdir,
            str(tmpdir.join('out.ff'))])

def test_build_constraints(tmpdir):
    subprocess.run(["viparr-build-constraints", f"{dmsdir}/ww.dms", "-o", f"{tmpdir}/out.dms"], check=True)

def test_compare_forcefields(tmpdir):
    subprocess.run(["viparr-compare-forcefields", f"{ffdir}/amber03", f"{ffdir}/amber99"], check=True)

