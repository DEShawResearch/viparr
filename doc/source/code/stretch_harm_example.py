import msys
import viparr

def stretch_harm(tsystem, ff):
    if len(ff.getParams('stretch_harm')) == 0:
        raise RuntimeError('Plugin "bonds" requires "stretch_harm" params')
    if len(tsystem.nonPseudoBonds) == 0:
        return
    matcher = viparr.ParameterMatcher.FromFF(ff, 'stretch_harm',
            viparr.SystemToPattern.Bonded,
            viparr.TypeToPattern.Default,
            [viparr.Permutation.Identity, viparr.Permutation.Reverse])
    table = tsystem.system.addTable('stretch_harm', 2,
            viparr.Forcefield.ParamTable('stretch_harm'))
    table.category = 'bond'
    table.addTermProp('constrained', int)
    for bond in tsystem.nonPseudoBonds:
        param, perm = matcher.match(tsystem, bond)
        if param.id == -1:
            raise RuntimeError('"stretch_harm" parameter not found')
        table.addTerm(bond, param)

viparr.Forcefield.DelPlugin('bonds')
viparr.Forcefield.AddPlugin('bonds', viparr.Forcefield.Plugin(stretch_harm))

mol = msys.Load('my_system.dms')
ff = viparr.ImportForcefield('my_ff')
viparr.ExecuteViparr(mol, [ff])
msys.SaveDMS(mol, 'out.dms')
