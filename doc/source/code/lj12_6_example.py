import msys
import viparr
from math import sqrt

lj12_6 = viparr.Rules.VDWFunc()
lj12_6.vdw_table_name = 'vdw_12_6'
lj12_6.param_names = ['sigma', 'epsilon']
lj12_6.pair_table_name = 'pair_12_6_es'
lj12_6.pair_param_names = ['aij', 'bij']
lj12_6.supported_rules = ['geometric', 'arithmetic/geometric']
viparr.Rules.DelVDWFunc('lj12_6_sig_epsilon')
viparr.Rules.AddVDWFunc('lj12_6_sig_epsilon', lj12_6)

def geometric(vdw_a, vdw_b, scale_factor):
    sij = sqrt(vdw_a[0] * vdw_b[0])
    eij = sqrt(vdw_a[1] * vdw_b[1])
    aij = sij ** 12 * eij * 4.0 * scale_factor
    bij = sij ** 6 * eij * 4.0 * scale_factor
    return [aij, bij]
viparr.Rules.DelVDWCombRule('geometric')
viparr.Rules.AddVDWCombRule('geometric', geometric)

mol = msys.Load('my_system.dms')
ff = viparr.ImportForcefield('my_ff')
viparr.ExecuteViparr(mol, [ff])
msys.SaveDMS(mol, 'out.dms')
