#include <boost/python.hpp>
#include "version.hxx"
#include "../src/add_system_tables.hxx"
#include "../src/postprocess/build_constraints.hxx"
#include "../src/postprocess/fix_masses.hxx"
#include "../src/execute_viparr.hxx"
#include "../src/execute_iviparr.hxx"
#include "../src/util/get_bonds_angles_dihedrals.hxx"
#include "../src/util/system_to_dot.hxx"
#include <msys/version.hxx>

using namespace boost::python;
using namespace desres;
using namespace desres::viparr;

namespace desres { namespace viparr {

    void export_compile_plugins();
    void export_exports();
    void export_forcefield();
    void export_imports();
    void export_merges();
    void export_parameter_matcher();
    void export_pattern();
    void export_permutation();
    void export_rules();
    void export_system_to_pattern();
    void export_templated_system();
    void export_type_to_pattern();
    void export_typers();

}}

namespace {

    void add_system_tables(msys::SystemPtr sys, const object& py_share_params, const object& py_name_map) {
        dict D(py_share_params);
        list L(D.keys());
        std::map<std::string, std::list<msys::Id> > share_params;
        for (unsigned i = 0; i < len(L); ++i) {
            list py_vals(D[L[i]]);
            std::list<msys::Id> vals;
            for (unsigned j = 0; j < len(py_vals); ++j)
                vals.push_back(extract<msys::Id>(py_vals[j]));
            std::string key = extract<std::string>(L[i]);
            share_params.insert(std::make_pair(key, vals));
        }
        if (py_name_map == object())
            AddSystemTables(sys, share_params);
        else {
            dict D(py_name_map);
            std::map<std::string, std::string> name_map;
            list L = D.keys();
            for (unsigned i = 0; i < len(L); ++i) {
                std::string key = extract<std::string>(L[i]);
                std::string value = extract<std::string>(D[L[i]]);
                name_map.insert(std::make_pair(key, value));
            }
            AddSystemTables(sys, share_params, name_map);
        }
    }

    void build_constraints(msys::SystemPtr sys, const object& py_atoms, bool keep, const object& py_exclude, bool verbose) {
        list L(py_atoms);
        msys::IdList atoms;
        for (unsigned i = 0; i < len(L); ++i)
            atoms.push_back(extract<msys::Id>(py_atoms[i]));
        L = list(py_exclude);
        std::set<std::string> exclude;
        for (unsigned i = 0; i < len(L); ++i) {
            std::string excluded = extract<std::string>(L[i]);
            exclude.insert(excluded);
        }
        BuildConstraints(sys, atoms, keep, exclude, true, verbose); // Always optimize vsite defs
    }

    void execute_viparr(msys::SystemPtr sys, const object& py_fflist, const object& py_selection, bool rename_atoms, bool rename_residues, 
                        bool with_constraints, bool fix_masses, bool fatal, bool compile_plugins, bool verbose, bool verbose_matching) {
        list L(py_fflist);
        std::vector<ForcefieldPtr> fflist;
        for (unsigned i = 0; i < len(L); ++i) {
            ForcefieldPtr ff = extract<ForcefieldPtr>(L[i]);
            fflist.push_back(ff);
        }
        L = list(py_selection);
        msys::IdList selection;
        for (unsigned i = 0; i < len(L); ++i)
            selection.push_back(extract<msys::Id>(L[i]));
        ExecuteViparr(sys, fflist, selection, rename_atoms, rename_residues, with_constraints, true, fix_masses, fatal, compile_plugins, verbose, verbose_matching); // Always optimize vsite defs
    }

    ForcefieldPtr execute_iviparr(msys::SystemPtr sys, object py_atoms, RulesPtr rules, bool templateonly) {
        list L(py_atoms);
        msys::IdList atoms;
        for (unsigned i = 0; i < len(L); ++i) {
            atoms.push_back(extract<msys::Id>(L[i]));
        }
        return ExecuteIviparr(sys, atoms, rules, templateonly);
    }

    void fix_masses(msys::SystemPtr sys, const object& py_atoms, bool verbose) {
        list L(py_atoms);
        msys::IdList atoms;
        for (unsigned i = 0; i < len(L); ++i)
            atoms.push_back(extract<msys::Id>(L[i]));
        FixMasses(sys, atoms, verbose);
    }

    list convert_vector_idlist(std::vector<msys::IdList> const& v) {
        list L;
        for (auto ids : v) {
            list s;
            for (auto id : ids) {
                s.append(id);
            }
            L.append(s);
        }
        return L;
    }
    dict get_bonds_angles_dihedrals(msys::SystemPtr mol) {
        std::vector<msys::IdList> non_pseudo_bonds;
        std::vector<msys::IdList> pseudo_bonds;
        std::vector<msys::IdList> angles;
        std::vector<msys::IdList> dihedrals;
        desres::viparr::GetBondsAnglesDihedrals(mol, mol->atoms(), non_pseudo_bonds, pseudo_bonds, angles, dihedrals);
        dict d;
        d["non_pseudo_bonds"] = convert_vector_idlist(non_pseudo_bonds);
        d["pseudo_bonds"] = convert_vector_idlist(pseudo_bonds);
        d["angles"] = convert_vector_idlist(angles);
        d["dihedrals"] = convert_vector_idlist(dihedrals);
        return d;
    }

    std::string system_to_dot(msys::SystemPtr mol, msys::Id residue_id) {
        std::stringstream ss;
        SystemToDot(mol, ss, residue_id);
        return ss.str();
    }

}


BOOST_PYTHON_MODULE(_viparr) {
    boost::python::scope().attr("version")=std::string(VIPARR_VERSION);
    boost::python::scope().attr("hexversion")=VIPARR_VERSION_HEX;
    boost::python::scope().attr("msys_version")=std::string(MSYS_VERSION);
    def("BuildConstraints", build_constraints);
    def("AddSystemTables", add_system_tables);
    def("ExecuteViparr", execute_viparr);
    def("ExecuteIviparr", execute_iviparr);
    def("FixMasses", fix_masses);
    def("ReorderIDs", desres::viparr::ReorderIDs);
    def("GetBondsAnglesDihedrals", get_bonds_angles_dihedrals);
    def("SystemToDot", system_to_dot);
    desres::viparr::export_compile_plugins();
    desres::viparr::export_exports();
    desres::viparr::export_forcefield();
    desres::viparr::export_imports();
    desres::viparr::export_merges();
    desres::viparr::export_parameter_matcher();
    desres::viparr::export_pattern();
    desres::viparr::export_permutation();
    desres::viparr::export_rules();
    desres::viparr::export_system_to_pattern();
    desres::viparr::export_templated_system();
    desres::viparr::export_type_to_pattern();
    desres::viparr::export_typers();
}
