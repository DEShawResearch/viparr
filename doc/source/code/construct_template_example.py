import msys
import viparr

# Load the system, create a clone with just the ligand, and create a template
mol = msys.Load('my_ligand.dms')
mol = mol.clone('resname LIG')
tpl_mol = viparr.TemplatedSystem(mol)

# Set the atom names, charges, and atom types
for atom in mol.atoms:
    atom.name = 'name_%d' % atom.id # Make sure that names are unique
    atom.charge = 0.0
    tpl_mol.setTypes(atom, 'bonded_type', 'nonbonded_type')

# Add some impropers and extra exclusions
tpl_mol.addImproper([mol.atom(0), mol.atom(1), mol.atom(2), mol.atom(3)])
tpl_mol.addExclusion([mol.atom(0), mol.atom(1)])

# Write the template to disk
# Template name is the name of the residue of the first atom in the template
viparr.ExportTemplates([tpl_mol], 'template.my_ligand')
