#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "version.hxx"
#include "../src/add_system_tables.hxx"
#include "../src/ff.hxx"
#include "../src/importexport/import_ff.hxx"
#include "../src/importexport/export_ff.hxx"
#include "../src/pattern.hxx"
#include "../src/postprocess/build_constraints.hxx"
#include "../src/postprocess/optimize_vsitedefs.hxx"
#include "../src/postprocess/fix_masses.hxx"
#include "../src/execute_viparr.hxx"
#include "../src/execute_iviparr.hxx"
#include "../src/util/get_bonds_angles_dihedrals.hxx"
#include "../src/util/system_to_dot.hxx"
#include <msys/version.hxx>
#include <msys/python/capsule.hxx>

using namespace pybind11;
using namespace desres::msys;
using namespace desres::msys::python;
using namespace desres::viparr;

namespace pybind11 { namespace detail {
    template <> struct type_caster<SystemPtr> {
    public:
        PYBIND11_TYPE_CASTER(SystemPtr, _("SystemPtr"));
        bool load(handle src, bool) { value = system_from_capsule(src.ptr()); return bool(value); }
        static handle cast(SystemPtr src, return_value_policy, handle) { return system_as_capsule(src); }
    };
    template <> struct type_caster<TermTablePtr> {
    public:
        PYBIND11_TYPE_CASTER(TermTablePtr, _("TermTablePtr"));
        bool load(handle src, bool) { value = termtable_from_capsule(src.ptr()); return bool(value); }
        static handle cast(TermTablePtr src, return_value_policy, handle) { return termtable_as_capsule(src); }
    };
    template <> struct type_caster<ParamTablePtr> {
    public:
        PYBIND11_TYPE_CASTER(ParamTablePtr, _("ParamTablePtr"));
        bool load(handle src, bool) { value = paramtable_from_capsule(src.ptr()); return bool(value); }
        static handle cast(ParamTablePtr src, return_value_policy, handle) { return paramtable_as_capsule(src); }
    };
}}

namespace {
    template <class Obj> bool eq(const Obj& self, const Obj& other) { return self==other; }
    template <class Obj> bool ne(const Obj& self, const Obj& other) { return self!=other; }
    template <class Obj> unsigned long hash(const Obj& obj) { return reinterpret_cast<unsigned long>(obj.get()); }
}


static void export_forcefield(module m) {
    class_<Forcefield, ForcefieldPtr>(m, "Forcefield")
        .def("__eq__", eq<ForcefieldPtr>)
        .def("__ne__", ne<ForcefieldPtr>)
        .def("__hash__", hash<ForcefieldPtr>)
        .def(init(&Forcefield::create))
        .def_static("ClearPlugins", []() { Forcefield::PluginRegistry().clear(); })
        .def_static("HasParamTable", &Forcefield::HasParamTable)
        .def_static("ParamTable", &Forcefield::ParamTable)
        .def_static("AddParamTable", &Forcefield::AddParamTable)
        .def_static("ClearParamTables", &Forcefield::ClearParamTables)
        .def_static("AllParamTables", &Forcefield::AllParamTables)
        .def("rules", &Forcefield::rules)
        .def("resetRules", &Forcefield::resetRules)
        .def("typer", &Forcefield::typer)
        .def("resetTyper", &Forcefield::resetTyper)
        .def("paramTables", &Forcefield::paramTables)
        .def("rowIDs", &Forcefield::rowIDs)
        .def("cmapTable", &Forcefield::cmapTable)
        .def("cmapTables", &Forcefield::cmapTables)
        ;
}

static void export_pattern(module m) {
    class_<Pattern>(m, "Pattern")
        ;
    class_<SystemToPattern, SystemToPatternPtr>(m, "SystemToPattern")
        .def_readwrite_static("NBType", &SystemToPattern::NBType)
        .def_readwrite_static("BType", &SystemToPattern::BType)
        .def_readwrite_static("Bonded", &SystemToPattern::Bonded)
        .def_readwrite_static("BondToFirst", &SystemToPattern::BondToFirst)
        .def_readwrite_static("PseudoBType", &SystemToPattern::PseudoBType)
        .def_readwrite_static("PseudoBondToFirst", &SystemToPattern::PseudoBondToFirst)
        .def_readwrite_static("PseudoBondToSecond", &SystemToPattern::PseudoBondToSecond)
        ;
    class_<TypeToPattern, TypeToPatternPtr>(m, "TypeToPattern")
        .def_readwrite_static("Default", &TypeToPattern::Default)
        .def_readwrite_static("Pseudo", &TypeToPattern::Pseudo)
        ;
    class_<Permutation, PermutationPtr>(m, "Permutation")
        .def_readwrite_static("Identity", &Permutation::Identity)
        .def_readwrite_static("Reverse", &Permutation::Reverse)
        //.def_readwrite_static("Improper", &Permutation::Improper)
        //.def_readwrite_static("Null", &Permutation::Null)
        ;

}

static void export_rules(module m) {
    auto cls = class_<Rules, RulesPtr>(m, "Rules")
        .def("__eq__", eq<RulesPtr>)
        .def("__ne__", ne<RulesPtr>)
        .def("__hash__", hash<RulesPtr>)
        .def("es_scale", &Rules::es_scale)
        .def("lj_scale", &Rules::lj_scale)
        .def_readwrite("info", &Rules::info)
        .def_readwrite("vdw_func", &Rules::vdw_func)
        .def_readwrite("vdw_comb_rule", &Rules::vdw_comb_rule)
        .def_readwrite("nbfix_identifier", &Rules::nbfix_identifier)
        .def_readwrite("plugins", &Rules::plugins)
        .def("setExclusions", &Rules::setExclusions)
        .def_property_readonly("exclusions", &Rules::exclusions)
        .def(init(&Rules::create))
        .def_static("ClearVdwFuncRegistry", []() { Rules::VDWFuncRegistry().clear(); })
        .def_static("ClearVdwCombRuleRegistry", []() { Rules::VDWCombRuleRegistry().clear(); })
        .def_property_static("VDWFuncRegistry", [](object) { return Rules::VDWFuncRegistry(); },
                                                [](object, std::map<std::string, Rules::VDWFunc> r) { Rules::VDWFuncRegistry() = r; })
        .def_property_static("VDWCombRuleRegistry", [](object) { return Rules::VDWCombRuleRegistry(); },
                                                    [](object, std::map<std::string, Rules::VDWCombRulePtr> r) { Rules::VDWCombRuleRegistry() = r; })
        ;

    class_<Rules::VDWFunc>(cls, "VDWFunc")
        .def(init<>())
        .def_readwrite("vdw_table_name", &Rules::VDWFunc::vdw_table_name)
        .def_readwrite("param_names", &Rules::VDWFunc::param_names)
        .def_readwrite("pair_table_name", &Rules::VDWFunc::pair_table_name)
        .def_readwrite("pair_param_names", &Rules::VDWFunc::pair_param_names)
        ;
    class_<Rules::VDWCombRule, Rules::VDWCombRulePtr>(cls, "VDWCombRule")
        ;
}

static void export_typers(module m) {
    class_<TemplateTyper, TemplateTyperPtr>(m, "TemplateTyper")
        .def("__eq__", eq<RulesPtr>)
        .def("__ne__", ne<RulesPtr>)
        .def("__hash__", hash<RulesPtr>)
        .def(init(&TemplateTyper::create))
        .def("templates", &TemplateTyper::templates)
        .def("addTemplate", &TemplateTyper::addTemplate)
        .def("matchFragment", [](TemplateTyper& typer, TemplatedSystemPtr sys, IdList const& atoms) {
            std::stringstream ss;
            std::vector<std::pair<TemplatedSystemPtr, IdList> > matches;
            bool matched = typer.matchFragment(sys, atoms, matches, ss);
            if (!matched) {
                matches.clear();
            } else {
                ss.clear();
            }
            return make_tuple(cast(matches), str(ss.str()));
        })
        ;
    class_<TemplatedSystem, TemplatedSystemPtr>(m, "TemplatedSystem")
        .def("__eq__", eq<TemplatedSystemPtr>)
        .def("__ne__", ne<TemplatedSystemPtr>)
        .def("__hash__", hash<TemplatedSystemPtr>)
        .def(init([]() { return TemplatedSystem::create(); }))
        .def(init([](SystemPtr m) { return TemplatedSystem::create(m); }))
        .def("system", &TemplatedSystem::system)
        .def("btype", &TemplatedSystem::btype)
        .def("nbtype", &TemplatedSystem::nbtype)
        .def("pset", &TemplatedSystem::pset)
        .def("pseudoBonds", &TemplatedSystem::pseudoBonds)
        .def("nonPseudoBonds", &TemplatedSystem::nonPseudoBonds)
        .def("angles", &TemplatedSystem::angles)
        .def("dihedrals", &TemplatedSystem::dihedrals)
        .def("exclusions", &TemplatedSystem::exclusions)
        .def("impropers", &TemplatedSystem::impropers)
        .def("cmaps", &TemplatedSystem::cmaps)
        ;

}

PYBIND11_MODULE(_viparr, m) {
    m.add_object("version", str(VIPARR_VERSION));
    m.add_object("hexversion", int_(VIPARR_VERSION_HEX));
    m.add_object("msys_version", str(MSYS_VERSION));
    m.def("BuildConstraints", BuildConstraints);
    m.def("OptimizeVsiteDefs", OptimizeVsiteDefs);
    m.def("AddSystemTables", AddSystemTables);
    m.def("ExecuteViparr", ExecuteViparr);
    m.def("ExecuteIviparr", ExecuteIviparr);
    m.def("FixMasses", FixMasses);
    m.def("ReorderIDs", ReorderIDs);
    m.def("ImportForcefield", ImportForcefield);
    m.def("ImportRules", ImportRules);
    m.def("ImportTemplates", ImportTemplates);
    m.def("ImportCmap", ImportCmap);
    m.def("ImportParams", ImportParams);
    m.def("ExportForcefield", ExportForcefield);
    m.def("ExportRules", ExportRules);
    m.def("ExportTemplates", ExportTemplates);
    m.def("ExportCmap", ExportCmap);
    m.def("ExportParams", ExportParams);
    m.def("GetBondsAnglesDihedrals", [](SystemPtr mol) {
        std::vector<IdList> non_pseudo_bonds;
        std::vector<IdList> pseudo_bonds;
        std::vector<IdList> angles;
        std::vector<::IdList> dihedrals;
        GetBondsAnglesDihedrals(mol, mol->atoms(), non_pseudo_bonds, pseudo_bonds, angles, dihedrals);
        dict d;
        d["non_pseudo_bonds"] = cast(non_pseudo_bonds);
        d["pseudo_bonds"] = cast(pseudo_bonds);
        d["angles"] = cast(angles);
        d["dihedrals"] = cast(dihedrals);
        return d;
        });
    m.def("SystemToDot", [](SystemPtr mol, Id residue_id) {
        std::stringstream ss;
        SystemToDot(mol, ss, residue_id);
        return ss.str();
        });


    export_forcefield(m);
    export_pattern(m);
    export_rules(m);
    export_typers(m);
}


