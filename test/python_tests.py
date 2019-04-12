
import os, sys, unittest, platform
from viparr import *

class TestMain(unittest.TestCase):

    def testTemplatedSystem(self):
        tsys = TemplatedSystem()
        # Check initial system state and graph hash
        self.assertTrue(tsys.system.atoms == [])
        self.assertTrue(tsys.hash == '0')
        res = tsys.system.addResidue()
        res.name = 'tpl'
        self.assertTrue(tsys.name == 'tpl')
        tsys.name = 'tpl_new'
        self.assertTrue(res.name == 'tpl_new')
        a0 = res.addAtom()
        a1 = res.addAtom()
        a2 = res.addAtom()
        a3 = res.addAtom()
        a4 = res.addAtom()
        a5 = res.addAtom()
        a6 = res.addAtom()
        a7 = res.addAtom()
        a0.atomic_number = 1
        a1.atomic_number = 1
        a2.atomic_number = 1
        a3.atomic_number = 1
        a4.atomic_number = 1
        a5.atomic_number = 1
        a6.atomic_number = 1
        a7.atomic_number = 0
        # Check that graph hash is updated when atoms are added
        g = msys.Graph(tsys.system.atoms)
        self.assertTrue(tsys.hash == g.hash())
        b0 = a0.addBond(a1)
        b1 = a1.addBond(a2)
        a2.addBond(a3)
        a1.addBond(a4)
        a7.addBond(a0)
        a5.addBond(a0)
        a6.addBond(a0)
        g = msys.Graph(tsys.system.atoms)
        # Check that graph hash and graph are updated when bonds are added
        self.assertTrue(tsys.hash == g.hash())
        self.assertTrue(len(tsys.graph.match(g)) == 7)
        tsys.setTypes(a0, 'c', 'cnb', '')
        tsys.setTypes(a7, 'c_virt', 'c_virt', 'pset0')
        tsys.addTypedAtom(a0)
        tsys.addTypedAtom(a1)
        tsys.addTypedAtom(a2)
        tsys.addTypedAtom(a3)
        tsys.addTypedAtom(a4)
        tsys.addTypedAtom(a5)
        tsys.addTypedAtom(a6)
        tsys.addTypedAtom(a7)
        tsys.addNonPseudoBond([a0,a1])
        tsys.addNonPseudoBond([a1,a2])
        tsys.addNonPseudoBond([a2,a3])
        tsys.addNonPseudoBond([a1,a4])
        tsys.addPseudoBond([a0,a7])
        tsys.addAngle([a0,a1,a2])
        tsys.addAngle([a1,a2,a3])
        tsys.addAngle([a0,a1,a4])
        tsys.addAngle([a2,a1,a4])
        tsys.addDihedral([a0,a1,a2,a3])
        tsys.addDihedral([a4,a1,a2,a3])
        tsys.addImproper([a1,a0,a2,a4])
        tsys.addExclusion([a0,a3])
        tsys.addExclusion([a4,a3])
        tsys.addCmap([a0,a1,a2,a3,a4,a5,a6,a7])
        tsys.addPseudoType('some_type', 3)
        tsys.addPseudoSites('some_type', [a7,a0,a1])
        # Check cloning functionality
        tsys_copy1 = tsys.clone()
        tsys_copy2 = tsys.clone([a0.id, a1.id, a2.id])
        tsys_copy3 = tsys.clone('index 0 1 2')
        self.assertTrue(tsys_copy1 != tsys)
        self.assertTrue(tsys_copy1.hash == tsys.hash)
        self.assertTrue(tsys_copy2.hash == tsys_copy3.hash)
        # Check types, psets, tuple lists
        self.assertTrue(tsys.btype(a0) == 'c')
        self.assertTrue(tsys.nbtype(a0) == 'cnb')
        self.assertTrue(tsys.pset(a0) == '')
        self.assertTrue(tsys.btype(a7) == 'c_virt')
        self.assertTrue(tsys.nbtype(a7) == 'c_virt')
        self.assertTrue(tsys.pset(a7) == 'pset0')
        self.assertTrue(tsys_copy1.btype(a0) == 'c')
        self.assertTrue(tsys_copy1.nbtype(a0) == 'cnb')
        self.assertTrue(tsys_copy1.pset(a0) == '')
        self.assertTrue(tsys_copy1.btype(a7) == 'c_virt')
        self.assertTrue(tsys_copy1.nbtype(a7) == 'c_virt')
        self.assertTrue(tsys_copy1.pset(a7) == 'pset0')
        self.assertTrue(tsys_copy2.btype(a0) == 'c')
        self.assertTrue(tsys_copy2.nbtype(a0) == 'cnb')
        self.assertTrue(tsys_copy2.pset(a0) == '')
        self.assertTrue(tsys_copy3.btype(a0) == 'c')
        self.assertTrue(tsys_copy3.nbtype(a0) == 'cnb')
        self.assertTrue(tsys_copy3.pset(a0) == '')
        self.assertTrue(tsys.typedAtoms == [[a0],[a1],[a2],[a3],[a4],[a5],[a6],[a7]])
        self.assertTrue(len(tsys_copy2.typedAtoms) == 3)
        self.assertTrue(len(tsys_copy1.nonPseudoBonds) == len(tsys.nonPseudoBonds))
        self.assertTrue(len(tsys_copy2.nonPseudoBonds) == 2)
        tsys.removeTypedAtom(a7)
        tsys.removeTypedAtom(a1)
        self.assertTrue(tsys.typedAtoms == [[a0],[a2],[a3],[a4],[a5],[a6]])
        self.assertTrue(len(tsys_copy1.typedAtoms) == 8)
        self.assertRaises(RuntimeError, tsys.removeTypedAtom, a1)
        self.assertTrue(tsys.nonPseudoBonds == [[a0,a1],[a1,a2],[a2,a3],[a1,a4]])
        tsys.removeNonPseudoBond([a1,a2])
        self.assertTrue(tsys.nonPseudoBonds == [[a0,a1],[a2,a3],[a1,a4]])
        self.assertTrue(len(tsys_copy1.nonPseudoBonds) == 4)
        self.assertRaises(RuntimeError, tsys.removeNonPseudoBond, [a0,a4])
        self.assertTrue(tsys.pseudoBonds == [[a0,a7]])
        tsys.removePseudoBond([a0,a7])
        self.assertTrue(tsys.pseudoBonds == [])
        self.assertTrue(len(tsys_copy1.pseudoBonds) == 1)
        self.assertRaises(RuntimeError, tsys.removePseudoBond, [a1,a2])
        self.assertTrue(tsys.angles == [[a0,a1,a2],[a1,a2,a3],[a0,a1,a4],[a2,a1,a4]])
        tsys.removeAngle([a0,a1,a2])
        self.assertTrue(tsys.angles == [[a1,a2,a3],[a0,a1,a4],[a2,a1,a4]])
        self.assertTrue(len(tsys_copy1.angles) == 4)
        self.assertRaises(RuntimeError, tsys.removeAngle, [a2,a3,a1])
        self.assertTrue(tsys.dihedrals == [[a0,a1,a2,a3],[a4,a1,a2,a3]])
        tsys.removeDihedral([a4,a1,a2,a3])
        self.assertTrue(tsys.dihedrals == [[a0,a1,a2,a3]])
        self.assertTrue(len(tsys_copy1.dihedrals) == 2)
        self.assertRaises(RuntimeError, tsys.removeDihedral, [a0,a2,a1,a3])
        self.assertTrue(tsys.impropers == [[a1,a0,a2,a4]])
        tsys.removeImproper([a1,a0,a2,a4])
        self.assertTrue(tsys.impropers == [])
        self.assertTrue(len(tsys_copy1.impropers) == 1)
        self.assertRaises(RuntimeError, tsys.removeImproper, [a2,a0,a1,a4])
        self.assertTrue(tsys.exclusions == [[a0,a3],[a4,a3]])
        tsys.removeExclusion([a0,a3])
        self.assertTrue(tsys.exclusions == [[a4,a3]])
        tsys.removeExclusion([a3,a4]) # argument in reversed order
        self.assertTrue(tsys.exclusions == [])
        self.assertTrue(len(tsys_copy1.exclusions) == 2)
        self.assertRaises(RuntimeError, tsys.removeExclusion, [a3,a0])
        self.assertTrue(tsys.cmaps == [[a0,a1,a2,a3,a4,a5,a6,a7]])
        tsys.removeCmap([a0,a1,a2,a3,a4,a5,a6,a7])
        self.assertTrue(tsys.cmaps == [])
        self.assertTrue(len(tsys_copy1.cmaps) == 1)
        self.assertRaises(RuntimeError, tsys.removeCmap, [a0,a1,a2,a3,a4,a5,a6,a7])
        pseudo_types = tsys.pseudoTypes
        self.assertTrue(len(pseudo_types) == 1)
        self.assertTrue(pseudo_types[0].name == 'some_type')
        self.assertTrue(pseudo_types[0].nsites == 3)
        self.assertTrue(pseudo_types[0].sites_list == [[a7,a0,a1]])
        a3.remove()
        # Check that hash and graph are updated when atom and bonds are deleted
        g = msys.Graph(tsys.system.atoms)
        self.assertTrue(tsys.hash == g.hash())
        self.assertTrue(len(tsys.graph.match(g)) == 6)
        sys = tsys.system
        tsys2 = TemplatedSystem(sys)
        # Check constructor
        self.assertTrue(tsys.system == tsys2.system)
        self.assertTrue(tsys != tsys2)

    def testTemplateTyper(self):
        typer = TemplateTyper()
        self.assertTrue(typer.templates == [])
        tpl = TemplatedSystem()
        res = tpl.system.addResidue()
        res.name = 'my_template'
        a0 = res.addAtom()
        a1 = res.addAtom()
        a2 = res.addAtom()
        a3 = res.addAtom()
        a4 = res.addAtom()
        a0.addBond(a1)
        a1.addBond(a2)
        a1.addBond(a3)
        a0.addBond(a4)
        a0.name = 'a0'
        a1.name = 'a1'
        a2.name = 'a2'
        a3.name = 'a3'
        a4.name = 'pseudo'
        a0.charge = -0.5
        a1.charge = 0
        a2.charge = 0.25
        a3.charge = 0.25
        a4.charge = 0
        a0.atomic_number = 8
        a1.atomic_number = 6
        a2.atomic_number = 1
        a3.atomic_number = 1
        a4.atomic_number = 0
        tpl.setTypes(a0, 'o', 'o_nb', '')
        tpl.setTypes(a1, 'c', 'c_nb', '')
        tpl.setTypes(a2, 'h', 'h_nb', '')
        tpl.setTypes(a3, 'h', 'h_nb', '')
        tpl.setTypes(a4, 'o_virt', 'o_virt', 'pset0')
        tpl.addTypedAtom(a0)
        tpl.addTypedAtom(a1)
        tpl.addTypedAtom(a2)
        tpl.addTypedAtom(a3)
        tpl.addTypedAtom(a4)
        tpl.addNonPseudoBond([a0,a1])
        tpl.addNonPseudoBond([a1,a2])
        tpl.addNonPseudoBond([a1,a3])
        tpl.addPseudoBond([a0,a4])
        tpl.addImproper([a1,a0,a2,a3])
        tpl.addExclusion([a2,a3])
        tpl.addPseudoSites('some_pseudo', [a4,a0])
        typer.addTemplate(tpl)
        # Test addTemplate
        self.assertTrue(typer.templates == [tpl])
        self.assertTrue(typer.findTemplate('my_template') == [tpl])
        self.assertTrue(typer.findTemplate('other_template') == [])
        sys = TemplatedSystem()
        res = sys.system.addResidue()
        b0 = res.addAtom()
        b1 = res.addAtom()
        b2 = res.addAtom()
        b3 = res.addAtom()
        b0.addBond(b1)
        b1.addBond(b2)
        b1.addBond(b3)
        b0.atomic_number = 8
        b1.atomic_number = 6
        b2.atomic_number = 1
        b3.atomic_number = 1
        # Test matchFragment and assignMatch with copied atom and residue names
        match, error_msg = typer.matchFragment(sys, [b0,b1,b2,b3])
        self.assertTrue(len(match) > 0)
        typer.assignMatch(sys, match, True, True)
        self.assertTrue(b0.charge == -0.5)
        self.assertTrue(b1.charge == 0)
        self.assertTrue(b2.charge == 0.25)
        self.assertTrue(b3.charge == 0.25)
        self.assertTrue(b0.name == 'a0')
        self.assertTrue(b1.name == 'a1')
        self.assertTrue(b2.name == 'a2')
        self.assertTrue(b3.name == 'a3')
        self.assertTrue(len(sys.system.atoms) == 5)
        self.assertTrue(b0.nbonds == 2)
        p = b0.bonded_atoms[0]
        if p == b1:
            p = b0.bonded_atoms[1]
        self.assertTrue(p.charge == 0)
        self.assertTrue(p.name == 'pseudo')
        self.assertTrue(sys.btype(b0) == 'o')
        self.assertTrue(sys.nbtype(b0) == 'o_nb')
        self.assertTrue(sys.btype(b1) == 'c')
        self.assertTrue(sys.nbtype(b1) == 'c_nb')
        self.assertTrue(sys.btype(b2) == 'h')
        self.assertTrue(sys.nbtype(b2) == 'h_nb')
        self.assertTrue(sys.btype(b3) == 'h')
        self.assertTrue(sys.nbtype(b3) == 'h_nb')
        self.assertTrue(sys.btype(p) == 'o_virt')
        self.assertTrue(sys.nbtype(p) == 'o_virt')
        self.assertTrue(sys.pset(p) == 'pset0')
        self.assertTrue([b0] in sys.typedAtoms)
        self.assertTrue([b1] in sys.typedAtoms)
        self.assertTrue([b2] in sys.typedAtoms)
        self.assertTrue([b3] in sys.typedAtoms)
        self.assertTrue([p] in sys.typedAtoms)
        self.assertTrue(len(sys.typedAtoms) == 5)
        self.assertTrue([b0,b1] in sys.nonPseudoBonds)
        self.assertTrue([b1,b2] in sys.nonPseudoBonds)
        self.assertTrue([b1,b3] in sys.nonPseudoBonds)
        self.assertTrue(len(sys.nonPseudoBonds) == 3)
        self.assertTrue(sys.pseudoBonds == [[b0,p]])
        self.assertTrue([b0,b1,b2] in sys.angles)
        self.assertTrue([b0,b1,b3] in sys.angles)
        self.assertTrue([b2,b1,b3] in sys.angles)
        self.assertTrue(len(sys.angles) == 3)
        self.assertTrue(len(sys.dihedrals) == 0)
        self.assertTrue(sys.impropers == [[b1,b0,b2,b3]])
        self.assertTrue(sys.exclusions == [[b2,b3]])
        self.assertTrue(len(sys.cmaps) == 0)
        self.assertTrue(len(sys.pseudoTypes) == 1)
        self.assertTrue(sys.pseudoTypes[0].name == 'some_pseudo')
        self.assertTrue(sys.pseudoTypes[0].nsites == 2)
        self.assertTrue(sys.pseudoTypes[0].sites_list == [[p,b0]])
        self.assertTrue(res.name == 'my_template')
        typer.delTemplate(tpl)
        self.assertTrue(typer.templates == [])

    def testRules(self):
        # Test VDWFuncRegistry
        registry = Rules.GetVDWFuncRegistry()
        self.assertTrue('lj12_6_sig_epsilon' in registry)
        self.assertTrue('exp_6x' in registry)
        self.assertTrue(registry['lj12_6_sig_epsilon'].vdw_table_name == 'vdw_12_6')
        self.assertTrue(registry['lj12_6_sig_epsilon'].param_names == ['sigma', 'epsilon'])
        self.assertTrue(registry['lj12_6_sig_epsilon'].pair_table_name == 'pair_12_6_es')
        self.assertTrue(registry['lj12_6_sig_epsilon'].pair_param_names == ['aij', 'bij'])
        self.assertTrue(registry['lj12_6_sig_epsilon'].supported_rules == ['geometric', 'arithmetic/geometric'])
        new_func = Rules.VDWFunc()
        new_func.vdw_table_name = 'vdw_new'
        new_func.param_names = ['a']
        new_func.pair_table_name = 'pair_new'
        new_func.pair_param_names = ['b', 'c']
        new_func.supported_rules = ['new_rule']
        Rules.AddVDWFunc('new_func', new_func)
        self.assertTrue('new_func' in Rules.GetVDWFuncRegistry())
        self.assertTrue(Rules.GetVDWFuncRegistry()['new_func'] == new_func)
        Rules.DelVDWFunc('new_func')
        self.assertTrue('new_func' not in Rules.GetVDWFuncRegistry())

        # Test VDWCombRuleRegistry
        registry = Rules.GetVDWCombRuleRegistry()
        self.assertTrue('arithmetic/geometric' in registry)
        self.assertTrue('geometric' in registry)
        self.assertTrue('lb/geometric' in registry)
        self.assertTrue(len(registry['geometric']([1,1],[2,2],0.5)) == 2)
        def new_rule(param_A, param_B, scaleself):
            param_comb = [param_A[0] + param_B[0], param_A[1] + param_B[1]]
            param_comb.append(scale)
            return param_comb
        Rules.AddVDWCombRule('new_rule', new_rule)
        self.assertTrue('new_rule' in Rules.GetVDWCombRuleRegistry())
        self.assertTrue(Rules.GetVDWCombRuleRegistry()['new_rule'] == new_rule)
        Rules.DelVDWCombRule('new_rule')
        self.assertTrue('new_rule' not in Rules.GetVDWCombRuleRegistry())

        # Test Rules object
        rules = Rules()
        self.assertTrue(rules.info == [])
        self.assertTrue(rules.vdw_func == '')
        self.assertTrue(rules.vdw_comb_rule == '')
        self.assertTrue(rules.plugins == [])
        self.assertTrue(rules.exclusions == 1)
        rules.info = ['my rules']
        rules.vdw_func = 'new_func'
        rules.vdw_comb_rule = 'new_rule'
        rules.plugins = ['A', 'B']
        rules.setExclusions(4, [0,0,0.5],[0,0.2,0.8])
        self.assertTrue(rules.info == ['my rules'])
        self.assertTrue(rules.vdw_func == 'new_func')
        self.assertTrue(rules.vdw_comb_rule == 'new_rule')
        self.assertTrue(rules.plugins == ['A', 'B'])
        self.assertTrue(rules.exclusions == 4)
        self.assertTrue(rules.es_scale(2) == 0)
        self.assertTrue(rules.es_scale(3) == 0)
        self.assertTrue(rules.es_scale(4) == 0.5)
        self.assertTrue(rules.lj_scale(2) == 0)
        self.assertTrue(rules.lj_scale(3) == 0.2)
        self.assertTrue(rules.lj_scale(4) == 0.8)

    def testPattern(self):
        p = Pattern()
        self.assertTrue(p.atoms == [])
        self.assertTrue(p.bonds == [])
        self.assertTrue(p.flags == [])
        p.atoms = ['c', 'h']
        p.bonds = ['-']
        p.flags = ['hello']
        p2 = Pattern()
        p2.atoms = ['c', 'h']
        p2.bonds = ['-']
        p2.flags = ['hello']
        self.assertTrue(p == p2)
        self.assertTrue(p.atoms == ['c', 'h'])
        self.assertTrue(p.bonds == ['-'])
        self.assertTrue(p.flags == ['hello'])
        p2.flags = ['hi']
        self.assertTrue(p != p2)

    def testSystemToPattern(self):
        tsys = TemplatedSystem()
        res = tsys.system.addResidue()
        a0 = res.addAtom()
        a1 = res.addAtom()
        a2 = res.addAtom()
        a3 = res.addAtom()
        a4 = res.addAtom()
        a0.addBond(a1)
        bond = a1.addBond(a2)
        bond2 = a1.addBond(a3)
        a0.addBond(a4)
        bond.order = 2
        bond2.order = 3
        tsys.setTypes(a0, 'a', 'a_nb', '')
        tsys.setTypes(a1, 'b', 'b_nb', '')
        tsys.setTypes(a2, 'c', 'c_nb', '')
        tsys.setTypes(a3, 'd', 'd_nb', '')
        tsys.setTypes(a4, 'p', 'p', 'pset')
        pat = SystemToPattern.NBType(tsys, [a0,a1,a2])
        self.assertTrue(pat.atoms == ['a_nb', 'b_nb', 'c_nb'])
        self.assertTrue(pat.bonds == [])
        self.assertTrue(pat.flags == [])
        pat = SystemToPattern.BType(tsys, [a0,a1,a2])
        self.assertTrue(pat.atoms == ['a', 'b', 'c'])
        self.assertTrue(pat.bonds == [])
        self.assertTrue(pat.flags == [])
        pat = SystemToPattern.Bonded(tsys, [a0,a1,a2])
        self.assertTrue(pat.atoms == ['a', 'b', 'c'])
        self.assertTrue(pat.bonds == ['-', '='])
        self.assertTrue(pat.flags == [])
        pat = SystemToPattern.BondToFirst(tsys, [a1,a0,a2,a3])
        self.assertTrue(pat.atoms == ['b', 'a', 'c', 'd'])
        self.assertTrue(pat.bonds == ['-', '=', '#'])
        self.assertTrue(pat.flags == [])
        pat = SystemToPattern.PseudoBType(tsys, [a4, a0, a1])
        self.assertTrue(pat.atoms == ['a', 'b'])
        self.assertTrue(pat.bonds == [])
        self.assertTrue(pat.flags == ['pset'])
        pat = SystemToPattern.PseudoBondToFirst(tsys, [a4, a0, a1])
        self.assertTrue(pat.atoms == ['a', 'b'])
        self.assertTrue(pat.bonds == ['-'])
        self.assertTrue(pat.flags == ['pset'])

    def testTypeToPattern(self):
        type = 'a - b ~ c = dstuff-stuff # e'
        pat = TypeToPattern.Default(type)
        self.assertTrue(pat.atoms == ['a', 'b', 'c', 'dstuff-stuff', 'e'])
        self.assertTrue(pat.bonds == ['-', '~', '=', '#'])
        self.assertTrue(pat.flags == [])
        type += ' f'
        pat = TypeToPattern.Pseudo(type)
        self.assertTrue(pat.atoms == ['a', 'b', 'c', 'dstuff-stuff', 'e'])
        self.assertTrue(pat.bonds == ['-', '~', '=', '#'])
        self.assertTrue(pat.flags == ['f'])

    def testPermutation(self):
        pat = Pattern()
        pat.atoms = ['a', 'b', 'c', 'd']
        pat.bonds = ['-', '=', '#']
        pat2 = Permutation.Identity(pat)
        self.assertTrue(pat2.atoms == ['a', 'b', 'c', 'd'])
        self.assertTrue(pat2.bonds == ['-', '=', '#'])
        pat2 = Permutation.Reverse(pat)
        self.assertTrue(pat2.atoms == ['d', 'c', 'b', 'a'])
        self.assertTrue(pat2.bonds == ['#', '=', '-'])
        pat2 = Permutation.Improper1(pat)
        self.assertTrue(pat2.atoms == ['a', 'b', 'd', 'c'])
        self.assertTrue(pat2.bonds == ['-', '#', '='])
        pat2 = Permutation.Improper2(pat)
        self.assertTrue(pat2.atoms == ['a', 'c', 'b', 'd'])
        self.assertTrue(pat2.bonds == ['=', '-', '#'])
        pat2 = Permutation.Improper3(pat)
        self.assertTrue(pat2.atoms == ['a', 'c', 'd', 'b'])
        self.assertTrue(pat2.bonds == ['=', '#', '-'])
        pat2 = Permutation.Improper4(pat)
        self.assertTrue(pat2.atoms == ['a', 'd', 'b', 'c'])
        self.assertTrue(pat2.bonds == ['#', '-', '='])
        pat2 = Permutation.Improper5(pat)
        self.assertTrue(pat2.atoms == ['a', 'd', 'c', 'b'])
        self.assertTrue(pat2.bonds == ['#', '=', '-'])

    def testSimpleParameterMatching(self):
        table = msys.CreateParamTable()
        table.addProp('type', str)
        p0 = table.addParam()
        p1 = table.addParam()
        p2 = table.addParam()
        p3 = table.addParam()
        p4 = table.addParam()
        p5 = table.addParam()
        p0['type'] = '* c'
        p1['type'] = 'c c'
        p2['type'] = '? car'
        p3['type'] = 'c car'
        p4['type'] = 'c car'
        p5['type'] = '* car'
        matcher = ParameterMatcher([p0,p1,p2,p3,p4,p5], SystemToPattern.Bonded, TypeToPattern.Default, [Permutation.Identity, Permutation.Reverse])
        self.assertTrue(matcher.params == [p0,p1,p2,p3,p4,p5])
        self.assertTrue(matcher.sys_to_pattern == SystemToPattern.Bonded)
        self.assertTrue(matcher.type_to_pattern == TypeToPattern.Default)
        tsys = TemplatedSystem()
        res = tsys.system.addResidue()
        a0 = res.addAtom()
        a1 = res.addAtom()
        bond = a0.addBond(a1)
        tsys.setTypes(a0, 'car', 'car', '')
        tsys.setTypes(a1, 'c', 'c', '')
        self.assertRaises(RuntimeError, matcher.match, tsys, [a0,a1])
        match = matcher.match(tsys, [a0,a1], True)
        self.assertTrue(match[0] == p3)
        self.assertTrue(match[1] == Permutation.Reverse)
        term_table = tsys.system.addTable('pair', 2, table)
        matcher.writeMultiple(p3, [a0,a1], term_table)
        self.assertTrue(term_table.nterms == 2)
        self.assertTrue(term_table.term(0).param == p3)
        self.assertTrue(term_table.term(1).param == p4)
        tsys.setTypes(a1, 'car', 'car', '')
        match = matcher.match(tsys, [a1,a0], True)
        self.assertTrue(match[0] == p5)
        tsys.setTypes(a0, '?', '?', '')
        match = matcher.match(tsys, [a0,a1], True)
        self.assertTrue(match[0] == p2)
        self.assertTrue(match[1] == Permutation.Identity)
        tsys.setTypes(a1, 'c', 'c', '')
        match = matcher.match(tsys, [a1,a0], True)
        self.assertTrue(match[0] == p0)
        # Test custom SystemToPattern, TypeToPattern, and Permutation
        def stp(tsystem, atoms):
            p = Pattern()
            p.atoms = [tsystem.btype(atoms[0]), tsystem.btype(atoms[1])]
            return p
        def ttp(type):
            tokens = type.split()
            p = Pattern()
            p.atoms = [tokens[0], tokens[1]]
            return p
        def perm(p):
            q = Pattern()
            q.atoms = [p.atoms[1], p.atoms[0]]
            q.bonds = p.bonds
            q.flags = p.flags
            return q
        matcher = ParameterMatcher([p0,p1,p2,p3,p4,p5], stp, ttp, [perm])
        self.assertTrue(matcher.sys_to_pattern == stp)
        self.assertTrue(matcher.type_to_pattern == ttp)
        tsys.setTypes(a0, 'car', 'car', '')
        tsys.setTypes(a1, 'c', 'c', '')
        match = matcher.match(tsys, [a0,a1], True)
        self.assertTrue(match[0] == p3)
        self.assertTrue(match[1] == perm)
        tsys.setTypes(a0, 'h', 'h', '')
        tsys.setTypes(a0, 'h', 'h', '')
        match = matcher.match(tsys, [a0,a1])
        self.assertTrue(match[0].id == -1)
        self.assertTrue(match[1] == None)

    def testForcefield(self):
        # Test global param tables
        Forcefield.ClearParamTables()
        self.assertTrue(Forcefield.AllParamTables() == [])
        table = msys.CreateParamTable()
        Forcefield.AddParamTable('new_table', table)
        self.assertTrue(Forcefield.AllParamTables() == ['new_table'])
        self.assertTrue(Forcefield.HasParamTable('new_table') == True)
        self.assertTrue(Forcefield.ParamTable('new_table') == table)
        table.addProp('type', str)
        table.addProp('value', float)
        table.addProp('int_value', int)
        p1 = table.addParam()
        p1['type'] = 'A'
        p1['value'] = 1.0
        p1['int_value'] = 1
        p2 = table.addParam()
        p2['type'] = 'B'
        p2['value'] = 2.0
        p2['int_value'] = 2

        # Test plugin registry
        registry = Forcefield.GetPluginRegistry()
        self.assertTrue('angles' in registry)
        self.assertTrue('bonds' in registry)
        self.assertTrue('cmap' in registry)
        self.assertTrue('exclusions' in registry)
        self.assertTrue('exclusions_and_scaled_pairs' in registry)
        self.assertTrue('impropers' in registry)
        self.assertTrue('mass' in registry)
        self.assertTrue('mass2' in registry)
        self.assertTrue('pairs_lj_scaled_14' in registry)
        self.assertTrue('propers' in registry)
        self.assertTrue('propers_allowmissing' in registry)
        self.assertTrue('pseudopol_fermi' in registry)
        self.assertTrue('ureybradley' in registry)
        self.assertTrue('vdw1' in registry)
        self.assertTrue('virtuals' in registry)
        self.assertTrue('virtuals_regular' in registry)
        tsys = TemplatedSystem()
        ff = Forcefield(Rules(), TemplateTyper())
        registry['impropers'].apply(tsys, ff)

        # Test user-defined plugin
        def my_plugin_func(tsys, ff):
            tsys.system.addTable('my_table', 2)
            return
        my_plugin = Forcefield.Plugin(my_plugin_func)
        self.assertTrue(my_plugin.prerequisites == [])
        my_plugin.prerequisites = ['angles']
        self.assertTrue(my_plugin.prerequisites == ['angles'])
        Forcefield.AddPlugin('my_plugin', my_plugin)
        self.assertTrue('my_plugin' in Forcefield.GetPluginRegistry())
        Forcefield.GetPluginRegistry()['my_plugin'].apply(tsys, ff)
        self.assertTrue('my_table' in tsys.system.table_names)
        Forcefield.DelPlugin('my_plugin')
        self.assertTrue('my_plugin' not in Forcefield.GetPluginRegistry())

        # Test Forcefield object
        rules = Rules()
        typer = TemplateTyper()
        ff = Forcefield(rules, typer)
        self.assertTrue(ff.name == '')
        ff.name = 'new_path'
        self.assertTrue(ff.name == 'new_path')
        self.assertTrue(ff.rules == rules)
        self.assertTrue(isinstance(ff.typer, TemplateTyper))
        self.assertTrue(ff.typer == typer)
        self.assertTrue(ff.params('new_table') == [])
        self.assertTrue(ff.paramTables == [])
        ff.appendParam('new_table', p1)
        ff.appendParam('new_table', p2)
        self.assertTrue(ff.paramTables == ['new_table'])
        self.assertTrue(ff.params('new_table') == [p1,p2])
        p3 = ff.appendParam('new_table', type='C', value=3.0, int_value=3)
        self.assertTrue(ff.params('new_table') == [p1, p2, p3])
        self.assertTrue(ff.findParams('new_table', type='C') == [p3])
        self.assertTrue(ff.findParams('new_table', value=2.0) == [p2])
        self.assertTrue(ff.findParams('new_table', int_value=1) == [p1])
        self.assertTrue(ff.findParams('new_table', int_value=1, type='C') == [])
        self.assertTrue(ff.findParams('new_table', int_value=1, type='A', value=1.0) == [p1])

        # Test ParameterMatcher.FromFF
        matcher = ParameterMatcher.FromFF(ff, 'new_table', SystemToPattern.BType, TypeToPattern.Default, [Permutation.Identity])
        self.assertTrue(matcher.params == [p1,p2,p3])

        # Test Forcefield object
        ff2 = ff.copy()
        self.assertTrue(ff2 != ff)
        self.assertTrue(ff.__str__() == ff2.__str__())
        self.assertTrue(ff.params('new_table') == ff2.params('new_table'))
        self.assertTrue(ff.rules == ff2.rules)
        self.assertTrue(ff.typer == ff2.typer)
        ff.delParams('new_table', p1)
        self.assertTrue(ff.params('new_table') == [p2,p3])
        self.assertTrue(ff2.params('new_table') == [p1,p2,p3])
        ff.replaceParam('new_table', p2, p1)
        self.assertTrue(ff.params('new_table') == [p1,p3])
        p4 = ff.replaceParam('new_table', p1, type='D', int_value=4, value=4.0)
        self.assertTrue(ff.params('new_table') == [p4,p3])
        ff.appendParam('new_table', p2)
        self.assertTrue(ff.params('new_table') == [p4,p3,p2])
        ff.delParams('new_table', int_value=3)
        self.assertTrue(ff.params('new_table') == [p4,p2])
        ff.delParams('new_table', type='D', int_value=2)
        self.assertTrue(ff.params('new_table') == [p4,p2])
        ff.delParams('new_table', type='D', int_value=4)
        self.assertTrue(ff.params('new_table') == [p2])
        ff.delParams('new_table', [p1,p2])
        self.assertTrue(ff.paramTables == [])
        self.assertTrue(ff.params('new_table') == [])

        # Test global param tables
        Forcefield.ClearParamTables()
        self.assertTrue(Forcefield.AllParamTables() == [])
        self.assertTrue(Forcefield.HasParamTable('new_table') == False)
        Forcefield.AddParamTable('new_table', table)
        Forcefield.ClearParamTables()
        self.assertTrue(Forcefield.AllParamTables() == [])

        # Test Forcefield object
        self.assertTrue(ff.cmap_tables == [])
        cmap = msys.CreateParamTable()
        ff.addCmapTable(cmap)
        self.assertTrue(ff.cmap_tables == [cmap])

    def testBuildConstraints(self):
        Forcefield.ClearParamTables()
        sys = msys.CreateSystem()
        res = sys.addResidue()
        a0 = res.addAtom()
        a1 = res.addAtom()
        a2 = res.addAtom()
        a3 = res.addAtom()
        a4 = res.addAtom()
        a5 = res.addAtom()
        a6 = res.addAtom()
        a0.addBond(a1)
        a0.addBond(a2)
        a3.addBond(a4)
        a3.addBond(a5)
        a3.addBond(a6)
        a0.atomic_number = 8
        a1.atomic_number = 1
        a2.atomic_number = 1
        a3.atomic_number = 7
        a4.atomic_number = 1
        a5.atomic_number = 1
        a6.atomic_number = 1
        sh_params = msys.CreateParamTable()
        Forcefield.AddParamTable('stretch_harm', sh_params)
        sh_params.addProp('r0', float)
        sh_params.addProp('fc', float)
        sh_params.addProp('type', str)
        sh_p = sh_params.addParam()
        ah_params = msys.CreateParamTable()
        Forcefield.AddParamTable('angle_harm', ah_params)
        ah_params.addProp('theta0', float)
        ah_params.addProp('fc', float)
        ah_params.addProp('type', str)
        ah_p = ah_params.addParam()
        sh_terms = sys.addTable('stretch_harm', 2, sh_params)
        sh_terms.addTermProp('constrained', int)
        sh_terms.addTerm([a0,a1], sh_p)
        sh_terms.addTerm([a0,a2], sh_p)
        sh_terms.addTerm([a3,a4], sh_p)
        sh_terms.addTerm([a3,a5], sh_p)
        sh_terms.addTerm([a3,a6], sh_p)
        ah_terms = sys.addTable('angle_harm', 3, ah_params)
        ah_terms.addTermProp('constrained', int)
        ah_terms.addTerm([a1,a0,a2], ah_p)
        ah_terms.addTerm([a4,a3,a5], ah_p)
        ah_terms.addTerm([a4,a3,a6], ah_p)
        ah_terms.addTerm([a5,a3,a6], ah_p)
        sys2 = msys.CloneSystem(sys.atoms)
        sys3 = msys.CloneSystem(sys.atoms)
        # Test build all constraints
        BuildConstraints(sys, verbose=False)
        self.assertTrue('constraint_hoh' in sys.table_names)
        self.assertTrue('constraint_ah3' in sys.table_names)
        self.assertTrue(sys.table('constraint_hoh').nterms == 1)
        self.assertTrue(sys.table('constraint_hoh').term(0).atoms == [a0, a1, a2])
        self.assertTrue(sys.table('constraint_ah3').nterms == 1)
        self.assertTrue(sys.table('constraint_ah3').term(0).atoms == [a3, a4, a5, a6])
        for term in sh_terms.terms:
            self.assertTrue(term['constrained'] == 1)
        self.assertTrue(ah_terms.term(0)['constrained'] == 1)
        self.assertTrue(ah_terms.term(1)['constrained'] == 0)
        self.assertTrue(ah_terms.term(2)['constrained'] == 0)
        self.assertTrue(ah_terms.term(3)['constrained'] == 0)
        # Test keep
        BuildConstraints(sys2, keep=True, verbose=False)
        self.assertTrue('constraint_hoh' in sys2.table_names)
        self.assertTrue('constraint_ah3' in sys2.table_names)
        sh_terms = sys2.table('stretch_harm')
        for term in sh_terms.terms:
            self.assertTrue(term['constrained'] == 0)
        ah_terms = sys2.table('angle_harm')
        for term in ah_terms.terms:
            self.assertTrue(term['constrained'] == 0)
        # Test build subset of constraints
        BuildConstraints(sys3, keep=True, exclude=['hoh'], verbose=False)
        self.assertTrue('constraint_hoh' not in sys3.table_names)
        self.assertTrue('constraint_ah3' in sys3.table_names)

    def testAddSystemTables(self):
        Forcefield.ClearParamTables()
        stretch_harm = msys.CreateParamTable()
        stretch_harm.addProp('type', str)
        stretch_harm.addProp('my_float', float)
        stretch_harm.addProp('my_int', int)
        stretch_harm.addProp('r0', float)
        stretch_harm.addProp('fc', float)
        stretch_harm.addProp('memo', str)
        p = stretch_harm.addParam()
        Forcefield.AddParamTable('stretch_harm', stretch_harm)
        sys = msys.LoadDMS('test/dms/ww.dms', False)
        sys.table('nonbonded').override_params.addProp('nbfix_identifier', str)
        or_p = sys.table('nonbonded').override_params.addParam()
        p0 = sys.table('nonbonded').params.param(0)
        p1 = sys.table('nonbonded').params.param(1)
        sys.table('nonbonded').setOverride(p0, p1, or_p)
        AddSystemTables(sys, {'angle_harm': 'angles'})
        self.assertTrue(Forcefield.HasParamTable('stretch_harm'))
        self.assertTrue(len(Forcefield.ParamTable('stretch_harm').params) > 1)
        self.assertTrue(Forcefield.ParamTable('stretch_harm').params[0] == p)
        self.assertTrue(Forcefield.HasParamTable('angles'))
        self.assertTrue(Forcefield.HasParamTable('nonbonded'))
        self.assertTrue(Forcefield.HasParamTable('constraint_hoh'))
        self.assertTrue(Forcefield.HasParamTable('vdw2'))
        self.assertTrue(sys.table('nonbonded').override_params == Forcefield.ParamTable('vdw2'))
        self.assertTrue(len(Forcefield.ParamTable('vdw2').params) == 1)

    def testCompilePlugins(self):
        # We smoke test the plugin compilations here; test for correctness
        # should be done by validation against viparr3
        # Test ureybradley_harm, mass, pairs_es_scaled, and pairs_lj_scaled
        Forcefield.ClearParamTables()
        charmm = ImportForcefield('test/ff3/charmm27')
        tip4p = ImportForcefield('test/ff3/tip4p')
        sys = msys.LoadDMS('test/dms/ww.dms', True)
        for a in sys.atoms:
            a.mass = 0
            a.charge = 0
        ExecuteViparr(sys, [charmm, tip4p], compile_plugins=False, verbose=False)
        self.assertTrue('ureybradley_harm' in sys.table_names)
        self.assertTrue('mass' in sys.table_names)
        self.assertTrue('pairs_12_6_es' not in sys.table_names)
        self.assertTrue(sys.atom(0).mass == 0)
        nstretch = sys.table('stretch_harm').nterms
        nureybradley = sys.table('ureybradley_harm').nterms
        CompilePlugins.CompileUreyBradley(sys)
        self.assertTrue(sys.table('stretch_harm').nterms == nstretch + nureybradley)
        CompilePlugins.CompileMasses(sys)
        self.assertTrue(sys.atom(0).mass != 0)
        pairs_table = CompilePlugins.AddPairsTable(sys)
        CompilePlugins.CompilePairsESScaled(sys, pairs_table)
        self.assertTrue('pair_12_6_es' in sys.table_names)
        self.assertTrue(sys.table('pair_12_6_es').term(0).param['qij'] != 0)
        self.assertTrue(sys.table('pair_12_6_es').term(0).param['aij'] == 0)
        self.assertTrue(sys.table('pair_12_6_es').term(0).param['bij'] == 0)
        CompilePlugins.CompilePairsLJScaled(sys, pairs_table)
        self.assertTrue(sys.table('pair_12_6_es').term(0).param['aij'] != 0)
        self.assertTrue(sys.table('pair_12_6_es').term(0).param['bij'] != 0)
        self.assertTrue(sys.nonbonded_info.vdw_funct == 'lj12_6_sig_epsilon')
        CompilePlugins.CleanupSystem(sys)
        self.assertTrue('ureybradley_harm' not in sys.table_names)
        self.assertTrue('mass' not in sys.table_names)
        self.assertTrue(sys.nonbonded_info.vdw_funct == 'vdw_12_6')
        sys = msys.LoadDMS('test/dms/ww.dms', True)
        for a in sys.atoms:
            a.mass = 0
            a.charge = 0
        ExecuteViparr(sys, [charmm, tip4p], compile_plugins=False, verbose=False)
        self.assertTrue('ureybradley_harm' in sys.table_names)
        self.assertTrue('mass' in sys.table_names)
        self.assertTrue('pairs_12_6_es' not in sys.table_names)
        self.assertTrue(sys.atom(0).mass == 0)
        nstretch = sys.table('stretch_harm').nterms
        nureybradley = sys.table('ureybradley_harm').nterms
        CompilePlugins.CompilePlugins(sys)
        self.assertTrue(sys.table('stretch_harm').nterms == nstretch + nureybradley)
        self.assertTrue(sys.atom(0).mass != 0)
        self.assertTrue('pair_12_6_es' in sys.table_names)
        self.assertTrue(sys.table('pair_12_6_es').term(0).param['qij'] != 0)
        self.assertTrue(sys.table('pair_12_6_es').term(0).param['aij'] != 0)
        self.assertTrue(sys.table('pair_12_6_es').term(0).param['bij'] != 0)
        self.assertTrue(sys.nonbonded_info.vdw_funct == 'lj12_6_sig_epsilon')
        CompilePlugins.CleanupSystem(sys)
        self.assertTrue('ureybradley_harm' not in sys.table_names)
        self.assertTrue('mass' not in sys.table_names)
        self.assertTrue(sys.nonbonded_info.vdw_funct == 'vdw_12_6')
        # Test improper_trig
        Forcefield.ClearParamTables()
        amber = ImportForcefield('test/ff3/amber03')
        tip3p = ImportForcefield('test/ff3/tip3p')
        sys = msys.LoadDMS('test/dms/ww.dms', True)
        for a in sys.atoms:
            a.mass = 0
            a.charge = 0
        ExecuteViparr(sys, [amber, tip3p], compile_plugins=False, verbose=False)
        self.assertTrue('improper_trig' in sys.table_names)
        nimproper = sys.table('improper_trig').nterms
        ndihedral = sys.table('dihedral_trig').nterms
        CompilePlugins.CompileImproperTrig(sys)
        self.assertTrue(sys.table('dihedral_trig').nterms == nimproper + ndihedral)
        CompilePlugins.CleanupSystem(sys)
        self.assertTrue('improper_trig' not in sys.table_names)
        # Test scaled_pair_overrides
        Forcefield.ClearParamTables()
        charmm27_ww = ImportForcefield('test/iviparr_charmm27_ww')
        tip4p = ImportForcefield('test/ff3/tip4p')
        sys = msys.LoadDMS('test/dms/ww.dms', True)
        for a in sys.atoms:
            a.mass = 0
            a.charge = 0
        ExecuteViparr(sys, [charmm27_ww, tip4p], compile_plugins=False, verbose=False)
        self.assertTrue('scaled_pair_overrides' in sys.table_names)
        # Took the following assertion out, as the pairs table is now added by execute viparr
        # self.assertTrue('pair_12_6_es' not in sys.table_names)
 
        pairs_table = CompilePlugins.AddPairsTable(sys)
        CompilePlugins.CompileScaledPairOverrides(sys, pairs_table)
        self.assertTrue('pair_12_6_es' in sys.table_names)
        CompilePlugins.CleanupSystem(sys)
        self.assertTrue('scaled_pair_overrides' not in sys.table_names)

    def testImportExport(self):
        import shutil
        import os
        if os.path.isdir('charmm_copy'):
            shutil.rmtree('charmm_copy')
        Forcefield.ClearParamTables()
        ff = ImportForcefield('test/ff3/charmm27')
        self.assertTrue(Forcefield.ParamTable('stretch_harm').nparams > 0)
        table_print = PrintParams(Forcefield.ParamTable('angle_harm').params)
        ExportForcefield(ff, 'charmm_copy')
        Forcefield.ClearParamTables()
        ff = ImportForcefield('charmm_copy')
        # Test that export and reimport preserves tables, rules, etc.
        self.assertTrue(PrintParams(Forcefield.ParamTable('angle_harm').params) == table_print)
        Forcefield.ClearParamTables()
        rules = ImportRules('test/ff3/charmm27/rules')
        ExportRules(rules, 'charmm_copy/rules2')
        rules2 = ImportRules('charmm_copy/rules2')
        self.assertTrue(rules.__str__() == rules2.__str__())
        templates = ImportTemplates('test/ff3/charmm27/templates')
        ExportTemplates(templates, 'charmm_copy/templates2')
        templates2 = ImportTemplates('charmm_copy/templates2')
        self.assertTrue(len(templates) == len(templates2))
        for tpl, tpl2 in zip(templates, templates2):
            self.assertTrue(tpl.hash == tpl2.hash)
        cmap = ImportCmap('test/ff3/charmm27/cmap')
        ExportCmap(cmap, 'charmm_copy/cmap2')
        cmap2 = ImportCmap('charmm_copy/cmap2')
        self.assertTrue(len(cmap) == len(cmap2))
        for c1, c2 in zip(cmap, cmap2):
            self.assertTrue(PrintParams(c1.params) == PrintParams(c2.params))
        params = ImportParams('stretch_harm', 'test/ff3/charmm27/stretch_harm')
        ExportParams(params, 'charmm_copy/stretch_harm2')
        params2 = ImportParams('stretch_harm', 'charmm_copy/stretch_harm2')
        self.assertTrue(PrintParams(params) == PrintParams(params2))
        Forcefield.ClearParamTables()
        shutil.rmtree('charmm_copy')

    def testMerge(self):
        Forcefield.ClearParamTables()
        src = ImportForcefield('test/ff3/charmm27')
        vdw_func = src.rules.vdw_func
        vdw_comb_rule = src.rules.vdw_comb_rule
        # Test merge rules
        rules = Rules()
        plugins = src.rules.plugins
        info = src.rules.info
        rules.plugins = ['my_plugin', 'bonds']
        rules.info = ['my_mod']
        MergeRules(src.rules, rules, verbose=False)
        self.assertTrue(src.rules.plugins == plugins + ['my_plugin'])
        self.assertTrue(src.rules.info == info + ['my_mod'])
        self.assertTrue(src.rules.vdw_func == vdw_func)
        self.assertTrue(src.rules.vdw_comb_rule == vdw_comb_rule)
        # Test merge templates
        my_tpl = TemplatedSystem()
        res = my_tpl.system.addResidue()
        a0 = res.addAtom()
        a0.atomic_number = 1
        res.name = 'my_tpl'
        typer = TemplateTyper()
        typer.addTemplate(my_tpl)
        my_tpl = TemplatedSystem()
        res = my_tpl.system.addResidue()
        a0 = res.addAtom()
        a0.atomic_number = 1
        res.name = 'ARG'
        typer.addTemplate(my_tpl)
        templates = src.typer.templates
        self.assertRaises(RuntimeError, MergeTemplates, src.typer, typer, append_only=True, verbose=False)
        MergeTemplates(src.typer, typer, verbose=False)
        self.assertTrue(len(src.typer.templates) == len(templates) + 1)
        self.assertTrue(len(src.typer.findTemplate('my_tpl')) == 1)
        self.assertTrue(len(src.typer.findTemplate('ARG')) == 1)
        self.assertTrue(src.typer.findTemplate('ARG')[0] == my_tpl)
        # Test merge params
        table = msys.CreateParamTable()
        table.addProp('type', str)
        p = table.addParam()
        Forcefield.AddParamTable('new_table', table)
        self.assertTrue(MergeParams([], [p], verbose=False) == [p])
        table = Forcefield.ParamTable('angle_harm')
        p1 = table.addParam()
        p1['type'] = 'test'
        p2 = table.addParam()
        p2['type'] = table.param(0)['type']
        p3 = src.params('angle_harm')[1]
        nparams = len(src.params('angle_harm'))
        self.assertRaises(RuntimeError, MergeParams, src.params('angle_harm'), [p2, p2], append_only=True, verbose=False)
        merged = MergeParams(src.params('angle_harm'), [p1, p2], verbose=False)
        self.assertTrue(len(merged) == nparams + 1)
        # New patch params come first
        self.assertTrue(merged[0] == p1)
        # Then come src params
        self.assertTrue(merged[1] == p2)
        self.assertTrue(merged[2] == p3)
        # Test merge forcefield with template typer
        patch = Forcefield(rules, typer)
        patch.appendParam('new_table', p)
        patch.appendParam('angle_harm', p1)
        patch.appendParam('angle_harm', p2)
        patch.addCmapTable(msys.CreateParamTable())
        MergeForcefields(src, patch, verbose=False)
        self.assertTrue(len(src.typer.templates) == len(templates) + 1)
        self.assertTrue(len(src.typer.findTemplate('my_tpl')) == 1)
        self.assertTrue(len(src.typer.findTemplate('ARG')) == 1)
        self.assertTrue(src.typer.findTemplate('ARG')[0] == my_tpl)
        self.assertTrue(src.params('new_table') == [p])
        self.assertTrue(len(src.params('angle_harm')) == nparams + 1)
        self.assertTrue(src.params('angle_harm')[0] == p1)
        self.assertTrue(src.params('angle_harm')[1] == p2)
        self.assertTrue(src.params('angle_harm')[2] == p3)
        self.assertTrue(len(src.cmap_tables) == 1)

    def testNBFix(self):
        Forcefield.ClearParamTables()
        lipids = ImportForcefield('test/ff3/charmm36_lipids')
        ions = ImportForcefield('test/ff3/charmm36_ions_nbfix')
        tip4p = ImportForcefield('test/ff3/tip4p')
        self.assertTrue(Forcefield.HasParamTable('vdw1'))
        self.assertTrue('nbfix_identifier' in Forcefield.ParamTable('vdw1').props)
        self.assertTrue(Forcefield.HasParamTable('vdw2'))
        self.assertTrue('nbfix_identifier' in Forcefield.ParamTable('vdw2').props)
        sys = msys.LoadDMS('test/dms/POPS_ions.dms', structure_only=True)
        ExecuteViparr(sys, [lipids, ions, tip4p], verbose=False)
        self.assertTrue('nonbonded' in sys.table_names)
        noverrides = len(sys.select('atomicnumber 11')) * len(sys.select('atomicnumber 17 or ((name O13A or name O13B) and resname POPS)'))
        self.assertEqual(sys.table('nonbonded').count_overrides(), noverrides)
        sys = msys.LoadDMS('test/dms/POPS_ions.dms', structure_only=True)
        ExecuteViparr(sys, [lipids, ions, tip4p], verbose=False, compile_plugins=False)
        
        # The following test was removed: execute_viparr now calls ApplyNBFix by default
        CompilePlugins.ApplyNBFix(sys)
        self.assertEqual(sys.table('nonbonded').count_overrides(), noverrides)
        # It's OK if multiple identical override rows match---just take one of them
        Forcefield.ParamTable('vdw2').params[-3].duplicate()
        sys = msys.LoadDMS('test/dms/POPS_ions.dms', structure_only=True)
        ExecuteViparr(sys, [lipids, ions, tip4p], verbose=False)
        self.assertEqual(sys.table('nonbonded').count_overrides(), noverrides)
        # If nbfix identifier of vdw1 row does not match the identifier of vdw2 row, override should not be included
        for p in Forcefield.ParamTable('vdw1').params:
            p['nbfix_identifier'] = ''
        sys = msys.LoadDMS('test/dms/POPS_ions.dms', structure_only=True)
        ExecuteViparr(sys, [lipids, ions, tip4p], verbose=False)
        self.assertEqual(sys.table('nonbonded').count_overrides(), 0)
        # It's an error if multiple override rows have same type but different params
        new_p = Forcefield.ParamTable('vdw2').params[-1].duplicate()
        Forcefield.ParamTable('vdw2').params[-1]['sigma'] = 0
        ions.appendParam('vdw2', new_p)
        sys = msys.LoadDMS('test/dms/POPS_ions.dms', structure_only=True)
        self.assertRaises(RuntimeError, ExecuteViparr, sys, [lipids, ions, tip4p], verbose=False)

    def testViparrExtensibility(self):
        Forcefield.ClearParamTables()
        amber99 = ImportForcefield('test/ff3/amber99')
        tip4p = ImportForcefield('test/ff3/tip4p')
        Forcefield.AddParamTable('my_table', msys.CreateParamTable())
        # Test user-defined plugin
        def my_plugin_func(tsys, ff):
            tsys.system.addTable('my_table', 2, Forcefield.ParamTable('my_table'))
            tsys.system.table('my_table').category = 'bond'
            tsys.system.table('my_table').addTerm([tsys.system.atom(0), tsys.system.atom(1)])
            return
        my_plugin = Forcefield.Plugin(my_plugin_func)
        Forcefield.AddPlugin('my_plugin', my_plugin)
        # Test GenerateStandardPlugin
        plugin_stretchharm = GenerateStandardPlugin('stretch_harm', 'bonds')
        Forcefield.AddPlugin('stretch_harm', plugin_stretchharm)
        # Test user-defined VDW function
        new_func = Rules.VDWFunc()
        new_func.vdw_table_name = 'vdw_new'
        new_func.param_names = ['sigma', 'epsilon']
        new_func.pair_table_name = 'pair_new'
        new_func.pair_param_names = ['sigma_sum', 'epsilon_sum', 'scale']
        new_func.supported_rules = ['new_rule']
        Rules.AddVDWFunc('new_func', new_func)
        # Test user-defined combine rule
        def new_rule(param_A, param_B, scale):
            param_comb = [param_A[0] + param_B[0], param_A[1] + param_B[1]]
            param_comb.append(scale)
            return param_comb
        Rules.AddVDWCombRule('new_rule', new_rule)
        amber99.rules.plugins += ['my_plugin']
        amber99.rules.plugins += ['stretch_harm']
        amber99.rules.vdw_func = 'new_func'
        amber99.rules.vdw_comb_rule = 'new_rule'
        tip4p.rules.plugins += ['stretch_harm']
        sys = msys.LoadDMS('test/dms/ww.dms', structure_only=True)
        bond_count = sys.nbonds
        ExecuteViparr(sys, [amber99, tip4p], with_constraints=False, verbose=False)
        sys = ReorderIDs(sys)
        self.assertTrue('my_table' in sys.table_names)
        self.assertTrue(sys.table('my_table').params == Forcefield.ParamTable('my_table'))
        self.assertTrue('stretch_harm' in sys.table_names)
        self.assertTrue(sys.table('stretch_harm').params == Forcefield.ParamTable('stretch_harm'))
        self.assertTrue(sys.table('stretch_harm').nterms == bond_count * 2)
        self.assertTrue('nonbonded' in sys.table_names)
        self.assertTrue('sigma' in sys.table('nonbonded').params.props)
        self.assertTrue('epsilon' in sys.table('nonbonded').params.props)
        self.assertTrue(sys.nonbonded_info.vdw_funct == 'vdw_new')
        self.assertTrue('pair_new' in sys.table_names)
        self.assertTrue('sigma_sum' in sys.table('pair_new').params.props)
        self.assertTrue('epsilon_sum' in sys.table('pair_new').params.props)
        self.assertTrue('scale' in sys.table('pair_new').params.props)
        self.assertTrue('qij' in sys.table('pair_new').params.props)
        self.assertTrue(sys.nonbonded_info.vdw_rule == 'new_rule')
        self.assertTrue('exclusion' in sys.table_names)

    def testViparrParametrizeFragment(self):
        Forcefield.ClearParamTables()
        amber99 = ImportForcefield('test/ff3/amber99')
        tip4p = ImportForcefield('test/ff3/tip4p')
        sys = msys.LoadDMS('test/dms/ww.dms', structure_only=True)
        water_atoms = sys.select('water')
        non_water_atoms = sys.select('not water')
        ExecuteViparr(sys, [amber99, tip4p], atoms=water_atoms, verbose=False)
        self.assertTrue('angle_harm' in sys.table_names)
        self.assertTrue(sys.table('angle_harm').nterms == len(water_atoms) / 3)
        self.assertTrue('constraint_hoh' in sys.table_names)
        self.assertTrue(Forcefield.HasParamTable('constraint_hoh'))
        self.assertTrue(sys.table('constraint_hoh').nterms == len(water_atoms) / 3)
        self.assertTrue(sys.natoms == len(non_water_atoms) + len(water_atoms) / 3 * 4)
        ExecuteViparr(sys, [amber99], atoms=non_water_atoms, verbose=False)
        self.assertTrue(sys.natoms == len(non_water_atoms) + len(water_atoms) / 3 * 4)
        self.assertTrue('nonbonded' in sys.table_names)
        self.assertTrue(sys.table('nonbonded').nterms == sys.natoms)
        ExecuteViparr(sys, [amber99, tip4p], verbose=False)
        self.assertTrue(sys.table('nonbonded').nterms == sys.natoms)
        self.assertTrue(sys.table('constraint_hoh').nterms == len(water_atoms) / 3)
        self.assertTrue('constraint_ah2' in sys.table_names)
        nconstraints = sys.table('constraint_ah2').nterms
        ExecuteViparr(sys, [amber99], atoms=non_water_atoms, with_constraints=False, verbose=False)
        BuildConstraints(sys, atoms=non_water_atoms, verbose=False)
        self.assertTrue(sys.table('nonbonded').nterms == sys.natoms)
        self.assertTrue(sys.table('constraint_hoh').nterms == len(water_atoms) / 3)
        self.assertTrue(sys.table('constraint_ah2').nterms == nconstraints)

    def testFixMasses(self):
        Forcefield.ClearParamTables()
        amber99 = ImportForcefield('test/ff3/amber99')
        tip4p = ImportForcefield('test/ff3/tip4p')
        sys = msys.LoadDMS('test/dms/ww_solv.dms', structure_only=True)
        ExecuteViparr(sys, [amber99, tip4p], fix_masses=False, verbose=False)
        self.assertTrue('constraint_ah1' in sys.table_names)
        atoms = sys.select('atomicnumber 8')
        import numpy
        masses = numpy.array([atom.mass for atom in atoms])
        median = numpy.median(masses)
        self.assertTrue(len(numpy.nonzero(masses == median)[0]) < len(masses))
        FixMasses(sys, verbose=False)
        masses = numpy.array([atom.mass for atom in atoms])
        self.assertTrue(len(numpy.nonzero(masses == median)[0]) == len(masses))
        ExecuteViparr(sys, [amber99], atoms=sys.select('protein'), fix_masses=False, verbose=False)
        masses = numpy.array([atom.mass for atom in atoms])
        self.assertTrue(len(numpy.nonzero(masses == median)[0]) < len(masses))
        ExecuteViparr(sys, [amber99, tip4p], verbose=False)
        masses = numpy.array([atom.mass for atom in atoms])
        self.assertTrue(len(numpy.nonzero(masses == median)[0]) == len(masses))

    def testViparrNeutralize(self):
        from viparr import neutralize
        Forcefield.ClearParamTables()
        sys = msys.LoadDMS('test/dms/ww_solv.dms', structure_only=True)
        sys = neutralize.Neutralize(sys, concentration=0.0, verbose=False)
        self.assertTrue(len(sys.select('atomicnumber 11')) == 1)
        self.assertTrue(len(sys.select('atomicnumber 17')) == 0)
        sys = neutralize.Neutralize(sys, concentration=2.0, verbose=False)
        self.assertTrue(len(sys.select('atomicnumber 11')) == 85)
        self.assertTrue(len(sys.select('atomicnumber 17')) == 84)
        sys = neutralize.Neutralize(sys, concentration=0.0, verbose=False)
        self.assertTrue(len(sys.select('atomicnumber 11')) == 1)
        self.assertTrue(len(sys.select('atomicnumber 17')) == 0)
        sys = msys.LoadDMS('test/dms/ww_solv.dms', structure_only=False)
        Forcefield.ClearParamTables()
        charmm27 = 'test/ff3/charmm27'
        sys = neutralize.Neutralize(sys, concentration=0.0, ffdir=charmm27, verbose=False)
        self.assertTrue(len(sys.select('atomicnumber 11')) == 1)
        self.assertTrue(len(sys.select('atomicnumber 17')) == 0)
        sys.nonbonded_info.vdw_rule = "arithmetic/geometric"
        sys = neutralize.Neutralize(sys, concentration=2.0, ffdir=charmm27, verbose=False)
        self.assertTrue(len(sys.select('atomicnumber 11')) == 85)
        self.assertTrue(len(sys.select('atomicnumber 17')) == 84)
        self.assertTrue('nonbonded' in sys.table_names)
        nonbonded_atoms = [term.atoms[0] for term in sys.table('nonbonded').terms]
        self.assertTrue(sys.select('atomicnumber 11')[40] in nonbonded_atoms)
        self.assertTrue(sys.select('atomicnumber 17')[40] in nonbonded_atoms)

    def testViparrSolvate(self):
        import os
        if os.path.exists('watbox.dms'):
            os.remove('watbox.dms')
        from viparr import solvate
        Forcefield.ClearParamTables()
        watbox = solvate.Solvate(msys.CreateSystem(), verbose=False)
        self.assertTrue(len(watbox.select('water')) == len(watbox.atoms))
        sys = msys.LoadDMS('test/dms/ww_solv.dms', structure_only=True)
        sys = sys.clone('not water')
        self.assertTrue(len(sys.select('water')) == 0)
        sys = solvate.Solvate(sys, verbose=False)
        self.assertTrue(len(sys.select('water')) > 0)
        sys = msys.LoadDMS('test/dms/ww_solv.dms', structure_only=False)
        sys = sys.clone('not water')
        self.assertTrue(len(sys.select('water')) == 0)
        tip4p = 'test/ff3/tip4p'
        sys = solvate.Solvate(sys, ffdir=tip4p, verbose=False)
        self.assertTrue(len(sys.select('water')) > 0)
        self.assertTrue('nonbonded' in sys.table_names)
        nonbonded_atoms = [term.atoms[0] for term in sys.table('nonbonded').terms]
        self.assertTrue(sys.select('water')[0] in nonbonded_atoms)
        watbox = solvate.Solvate(msys.CreateSystem(), ffdir=tip4p, verbose=False)
        self.assertTrue(len(watbox.select('water')) == len(watbox.atoms))
        msys.SaveDMS(watbox, 'watbox.dms')
        sys = msys.LoadDMS('test/dms/ww_solv.dms', structure_only=False)
        sys = sys.clone('not water')
        self.assertTrue(len(sys.select('water')) == 0)
        sys = solvate.Solvate(sys, 'watbox.dms', verbose=False)
        self.assertTrue(len(sys.select('water')) > 0)
        self.assertTrue('nonbonded' in sys.table_names)
        nonbonded_atoms = [term.atoms[0] for term in sys.table('nonbonded').terms]
        self.assertTrue(sys.select('water')[0] in nonbonded_atoms)
        os.remove('watbox.dms')

if __name__ == "__main__":
    unittest.main(verbosity=2)
