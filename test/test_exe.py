from __future__ import print_function

import os
dmsdir='%s/dms' % os.path.dirname(__file__)

def test_iviparr(tmpdir):
    import garden
    import subprocess
    garden.load('viparr-ff/2.1.0c7/data')
    subprocess.check_call(['iviparr', '--ffname', 'aa.amber.ff99SB', '%s/e.dms' % dmsdir,
            str(tmpdir.join('out.ff'))])


